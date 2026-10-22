//
// Namespace CoreFitSC
//   Core selection and far-field model for a space-charge mesh fitted to the bunch core
//   (OPAL-cycl, FFT solver).
//
// Created by: Daniel Winklehner, MIT (2026)
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#include "Algorithms/CoreFitSC.h"

#include "Message/GlobalComm.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

namespace {
    /// Sums over the live particles within halfWidth of centre in every axis:
    /// count, charge, sum (r - centre) (3), sum (r - centre)^2 (3).
    void accumulateClipped(const Vector_t* R, const double* Q, const int* bin,
                           std::size_t localNum, const Vector_t& centre,
                           const Vector_t& halfWidth, double sums[8]) {
        std::fill(sums, sums + 8, 0.0);
        for (std::size_t i = 0; i < localNum; ++i) {
            if (bin[i] < 0) {
                continue;
            }
            const Vector_t dr = R[i] - centre;
            if (std::abs(dr[0]) > halfWidth[0] ||
                std::abs(dr[1]) > halfWidth[1] ||
                std::abs(dr[2]) > halfWidth[2]) {
                continue;
            }
            sums[0] += 1.0;
            sums[1] += Q[i];
            for (unsigned int d = 0; d < 3; ++d) {
                sums[2 + d] += dr[d];
                sums[5 + d] += dr[d] * dr[d];
            }
        }
    }

    /// Mean and rms from the reduced sums of accumulateClipped() (count > 0)
    void meanAndSigma(const double sums[8], const Vector_t& centre,
                      Vector_t& mean, Vector_t& sigma) {
        for (unsigned int d = 0; d < 3; ++d) {
            const double shift = sums[2 + d] / sums[0];
            mean[d] = centre[d] + shift;
            sigma[d] = std::sqrt(std::max(0.0, sums[5 + d] / sums[0] - shift * shift));
        }
    }

    /// Half-widths k * sigma_d of a test, at least flatWidth
    Vector_t getHalfWidth(double k, const Vector_t& sigma) {
        Vector_t width;
        for (unsigned int d = 0; d < 3; ++d) {
            width[d] = std::max(k * sigma[d], CoreFitSC::flatWidth);
        }
        return width;
    }
}


void CoreFitSC::IpplReducer::sum(double* values, unsigned int count) const {
    allreduce(values, count, std::plus<double>());
}


void CoreFitSC::IpplReducer::min(double* values, unsigned int count) const {
    allreduce(values, count, std::less<double>());
}


double CoreFitSC::CoreSelection::getFarChargeFraction() const {
    return chargeLive != 0.0 ? chargeFar / chargeLive : 0.0;
}


CoreFitSC::CoreSelection CoreFitSC::selectCore(const Vector_t* R, const double* Q,
                                               const int* bin, std::size_t localNum,
                                               double nSigma, double clip,
                                               const Reducer& reducer) {
    CoreSelection selection;
    double sums[8];

    // Pass 0: the plain mean and rms of the live particles, in one sum about the origin.
    // The local space-charge frame is centred on the bunch. In a frame that is not, the
    // cancellation only blurs the first clip width, as long as the rms exceeds about 1e-8 of
    // the offset: every later pass is summed about the mean of the pass before. In an axis in
    // which the bunch is flat away from the origin, this rms is round-off: the first clip
    // still keeps the plane, and the stop test below does not depend on it.
    const Vector_t unlimited(std::numeric_limits<double>::infinity());
    accumulateClipped(R, Q, bin, localNum, Vector_t(0.0), unlimited, sums);
    reducer.sum(sums, 8);
    selection.numLive = static_cast<std::size_t>(sums[0]);
    selection.chargeLive = sums[1];

    Vector_t mean(0.0), sigma(0.0);
    if (sums[0] > 0.0) {
        meanAndSigma(sums, Vector_t(0.0), mean, sigma);

        // Passes 1 to maxClipPasses: the particles within clip * sigma of the previous
        // pass, summed about its mean.
        double numPrevious = sums[0];
        for (unsigned int pass = 1; pass <= maxClipPasses; ++pass) {
            accumulateClipped(R, Q, bin, localNum, mean, getHalfWidth(clip, sigma), sums);
            reducer.sum(sums, 8);
            if (sums[0] == 0.0) {
                break;
            }
            Vector_t newMean, newSigma;
            meanAndSigma(sums, mean, newMean, newSigma);

            // an axis whose new rms is at most flatWidth has settled, whatever its rms
            // before, which after pass 0 can be round-off above flatWidth
            bool sigmaSettled = true;
            for (unsigned int d = 0; d < 3; ++d) {
                if (newSigma[d] > flatWidth &&
                    !(std::abs(newSigma[d] - sigma[d]) < clipTolerance * newSigma[d])) {
                    sigmaSettled = false;
                }
            }
            mean = newMean;
            sigma = newSigma;
            selection.passes = pass;
            if (sums[0] == numPrevious || sigmaSettled) {
                selection.converged = true;
                break;
            }
            numPrevious = sums[0];
        }
    }
    selection.mean = mean;
    selection.sigma = sigma;

    // Final pass: the core within nSigma * sigma in every axis, its exact bounds, the far
    // particles and the halo proper.
    const Vector_t coreWidth = getHalfWidth(nSigma, sigma);
    const Vector_t haloWidth = getHalfWidth(haloSigma, sigma);
    double extrema[7];      // min x, y, z; -max x, y, z; -(largest far distance)
    std::fill(extrema, extrema + 7, std::numeric_limits<double>::max());
    double farSums[3] = {0.0, 0.0, 0.0};    // far count, far charge, halo count
    for (std::size_t i = 0; i < localNum; ++i) {
        if (bin[i] < 0) {
            continue;
        }
        const Vector_t dr = R[i] - mean;
        bool isCore = true;
        bool isHalo = false;
        for (unsigned int d = 0; d < 3; ++d) {
            const double dist = std::abs(dr[d]);
            if (dist > coreWidth[d]) {
                isCore = false;
            }
            if (dist > haloWidth[d]) {
                isHalo = true;
            }
        }
        if (isHalo) {
            farSums[2] += 1.0;
        }
        if (isCore) {
            for (unsigned int d = 0; d < 3; ++d) {
                extrema[d]     = std::min(extrema[d], R[i][d]);
                extrema[3 + d] = std::min(extrema[3 + d], -R[i][d]);
            }
            continue;
        }
        farSums[0] += 1.0;
        farSums[1] += Q[i];
        double distance = 0.0;
        for (unsigned int d = 0; d < 3; ++d) {
            if (sigma[d] > flatWidth) {
                distance = std::max(distance, std::abs(dr[d]) / sigma[d]);
            }
        }
        extrema[6] = std::min(extrema[6], -distance);
    }
    reducer.sum(farSums, 3);
    reducer.min(extrema, 7);

    selection.numFar = static_cast<std::size_t>(farSums[0]);
    selection.numCore = selection.numLive - selection.numFar;
    selection.numHalo = static_cast<std::size_t>(farSums[2]);
    selection.chargeFar = farSums[1];
    for (unsigned int d = 0; d < 3; ++d) {
        selection.boundsMin[d] = extrema[d];
        selection.boundsMax[d] = -extrema[3 + d];
    }
    selection.farMaxDistance = selection.numFar > 0 ? -extrema[6] : 0.0;

    return selection;
}


void CoreFitSC::getMeshBox(const Vector_t* R, const int* bin, std::size_t localNum,
                           const CoreSelection& selection, double nSigma,
                           const Reducer& reducer, Vector_t& boxMin, Vector_t& boxMax) {
    // the core test of the final pass of selectCore(), here for the particles flagged lost
    const Vector_t coreWidth = getHalfWidth(nSigma, selection.sigma);
    double extrema[6];      // min x, y, z; -max x, y, z
    for (unsigned int d = 0; d < 3; ++d) {
        extrema[d]     = selection.boundsMin[d];
        extrema[3 + d] = -selection.boundsMax[d];
    }
    for (std::size_t i = 0; i < localNum; ++i) {
        if (bin[i] >= 0) {
            continue;
        }
        const Vector_t dr = R[i] - selection.mean;
        if (std::abs(dr[0]) > coreWidth[0] ||
            std::abs(dr[1]) > coreWidth[1] ||
            std::abs(dr[2]) > coreWidth[2]) {
            continue;
        }
        for (unsigned int d = 0; d < 3; ++d) {
            extrema[d]     = std::min(extrema[d], R[i][d]);
            extrema[3 + d] = std::min(extrema[3 + d], -R[i][d]);
        }
    }
    reducer.min(extrema, 6);

    for (unsigned int d = 0; d < 3; ++d) {
        boxMin[d] = extrema[d];
        boxMax[d] = -extrema[3 + d];
    }
}


CoreFitSC::ParticleClass CoreFitSC::classify(int bin, const Vector_t& r,
                                             const Vector_t& boundsMin,
                                             const Vector_t& boundsMax,
                                             const Vector_t& meshMin,
                                             const Vector_t& meshMax) {
    if (isInside(r, boundsMin, boundsMax)) {
        return bin >= 0 ? ParticleClass::CORE : ParticleClass::INSIDE;
    }
    return isInside(r, meshMin, meshMax) ? ParticleClass::MESH : ParticleClass::FAR;
}


CoreFitSC::MomentSums::MomentSums(const Vector_t& centre, double gamma):
    centre_m(centre),
    gamma_m(gamma)
{
    sums_m.fill(0.0);
}


void CoreFitSC::MomentSums::add(double q, const Vector_t& r) {
    const Vector_t s = toRestFrame(r, centre_m, gamma_m);
    sums_m[0] += q;
    for (unsigned int i = 0; i < 3; ++i) {
        sums_m[1 + i] += q * s[i];
    }
    // q s_i s_j, i <= j: xx, xy, xz, yy, yz, zz
    unsigned int l = 4;
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = i; j < 3; ++j, ++l) {
            sums_m[l] += q * s[i] * s[j];
        }
    }
}


CoreFitSC::CoreMoments CoreFitSC::MomentSums::reduce(const Reducer& reducer) const {
    std::array<double, 10> sums = sums_m;
    reducer.sum(sums.data(), 10);

    CoreMoments moments;
    moments.charge = sums[0];
    moments.quadrupole.fill(Vector_t(0.0));
    if (sums[0] == 0.0) {
        return moments;
    }
    for (unsigned int i = 0; i < 3; ++i) {
        moments.centroid[i] = sums[1 + i] / sums[0];
    }
    // second moments about the centroid (xx, xy, xz, yy, yz, zz), then the traceless
    // quadrupole; each off-diagonal value is computed once, so the tensor is symmetric
    double M[6];
    for (unsigned int i = 0, l = 0; i < 3; ++i) {
        for (unsigned int j = i; j < 3; ++j, ++l) {
            M[l] = sums[4 + l] - sums[0] * moments.centroid[i] * moments.centroid[j];
        }
    }
    const double trace = M[0] + M[3] + M[5];
    for (unsigned int i = 0, l = 0; i < 3; ++i) {
        for (unsigned int j = i; j < 3; ++j, ++l) {
            const double Qij = 3.0 * M[l] - (i == j ? trace : 0.0);
            moments.quadrupole[i][j] = Qij;
            moments.quadrupole[j][i] = Qij;
        }
    }
    return moments;
}


void CoreFitSC::computeMultipoleField(const CoreMoments& moments,
                                      const std::vector<Vector_t>& position,
                                      std::size_t begin, std::size_t end, bool withQuadrupole,
                                      std::vector<Vector_t>& field,
                                      std::vector<double>& indicator) {
    const std::array<Vector_t, 3>& Q = moments.quadrupole;
    for (std::size_t i = begin; i < end; ++i) {
        const Vector_t d = position[i] - moments.centroid;
        const double r2 = dot(d, d);
        const double r3 = r2 * std::sqrt(r2);
        const Vector_t Qd({dot(Q[0], d), dot(Q[1], d), dot(Q[2], d)});
        const double dQd = dot(d, Qd);

        const Vector_t monopole = (moments.charge / r3) * d;
        const Vector_t quadrupole = (2.5 * dQd / (r2 * r2 * r3)) * d - Qd / (r2 * r3);

        const double monopoleNorm = euclidean_norm(monopole);
        indicator[i] = monopoleNorm > 0.0 ? euclidean_norm(quadrupole) / monopoleNorm : 0.0;
        field[i] = withQuadrupole ? monopole + quadrupole : monopole;
    }
}


void CoreFitSC::addDirectField(const std::vector<Vector_t>& target,
                               std::size_t begin, std::size_t end,
                               const std::vector<Vector_t>& source,
                               const std::vector<double>& charge,
                               double eps, std::vector<Vector_t>& field) {
    const double eps2 = eps * eps;
    const std::size_t numSources = source.size();
    for (std::size_t i = begin; i < end; ++i) {
        const Vector_t& t = target[i];
        Vector_t E(0.0);
        for (std::size_t j = 0; j < numSources; ++j) {
            const Vector_t d = t - source[j];
            const double d2 = dot(d, d);
            if (d2 == 0.0) {
                continue;
            }
            const double r2 = d2 + eps2;
            E += (charge[j] / (r2 * std::sqrt(r2))) * d;
        }
        field[i] += E;
    }
}


std::pair<std::size_t, std::size_t> CoreFitSC::getTableRange(std::size_t n, int rank,
                                                             int numRanks) {
    const std::size_t p = rank;
    const std::size_t numP = numRanks;
    return std::make_pair(p * n / numP, (p + 1) * n / numP);
}


std::size_t CoreFitSC::getSampleStride(std::size_t numTargets, std::size_t numSources,
                                       double maxPairs) {
    const double pairs = static_cast<double>(numTargets) * static_cast<double>(numSources);
    if (pairs <= maxPairs) {
        return 1;
    }
    return static_cast<std::size_t>(std::ceil(pairs / maxPairs));
}


bool CoreFitSC::isInCoreShell(const Vector_t& r, const Vector_t& lower, const Vector_t& upper) {
    for (unsigned int d = 0; d < 3; ++d) {
        const double centre = 0.5 * (lower[d] + upper[d]);
        const double halfWidth = 0.5 * (upper[d] - lower[d]);
        if (std::abs(r[d] - centre) > (1.0 - coreShellFraction) * halfWidth) {
            return true;
        }
    }
    return false;
}


CoreFitSC::ModeState CoreFitSC::decideMode(const ModeState& previous,
                                           std::size_t numCore, std::size_t numLive,
                                           double farChargeFraction, double maxFarFraction,
                                           unsigned int minDwell) {
    const double core = numCore;
    const double live = numLive;
    const unsigned int maxDwell = std::max(minDwell, 1u);
    // dwell 0: the run has not solved yet
    const unsigned int dwell = (previous.dwell == 0 ? maxDwell : previous.dwell);
    const bool dwellReached = (dwell >= minDwell);
    if (previous.mode == Mode::CORE) {
        // a broken-up core leaves at once, a large far charge only after the dwell
        if (core < std::max(minCoreParticles, minCoreFraction * live) ||
            (dwellReached && farChargeFraction > maxFarFraction)) {
            return ModeState{Mode::FULL, 1};
        }
    } else if (dwellReached &&
               farChargeFraction < returnFarFactor * maxFarFraction &&
               core > returnCoreFraction * live &&
               core >= minCoreParticles) {
        return ModeState{Mode::CORE, 1};
    }
    return ModeState{previous.mode, dwell < maxDwell ? dwell + 1 : maxDwell};
}
