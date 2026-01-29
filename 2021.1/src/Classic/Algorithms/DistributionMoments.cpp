//
// Class DistributionMoments
//   Computes the statistics of particle distributions.
//
// Copyright (c) 2021, Christof Metzger-Kraus
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "DistributionMoments.h"

#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include "Message/GlobalComm.h"
#include "Utility/Inform.h"

#include "OpalParticle.h"
#include "PartBunchBase.h"

extern Inform* gmsg;

DistributionMoments::DistributionMoments()
{
    reset();
}

void DistributionMoments::compute(PartBunchBase<double, 3> const& bunch)
{
    computeStatistics(bunch.begin(), bunch.end());
}

void DistributionMoments::compute(const std::vector<OpalParticle>::const_iterator &first,
                                  const std::vector<OpalParticle>::const_iterator &last)
{
    computeStatistics(first, last);
}

template<class InputIt>
void DistributionMoments::computeMeans(const InputIt &first, const InputIt &last)
{
    unsigned int localNum = last - first;
    std::vector<double> localStatistics(10);
    for (InputIt it = first; it != last; ++ it) {
        OpalParticle const& particle = *it;
        if (isParticleExcluded(particle)) {
            -- localNum;
            continue;
        }
        unsigned int l = 0;
        for (unsigned int i = 0; i < 6; ++ i, ++ l) {
            localStatistics[l] += particle[i];
        }

        // l = 6
        localStatistics[l++] += particle.getTime();

        double gamma = Util::getGamma(particle.getP());
        double eKin = (gamma - 1.0) * particle.getMass();
        localStatistics[l++] += eKin;
        localStatistics[l++] += gamma;
    }
    localStatistics.back() = localNum;

    allreduce(localStatistics.data(), localStatistics.size(), std::plus<double>());
    totalNumParticles_m = localStatistics.back();

    double perParticle = 1.0 / totalNumParticles_m;
    unsigned int l = 0;

    for (; l < 6; ++ l) {
        centroid_m[l] = localStatistics[l];
    }
    for (unsigned int i = 0; i < 3; ++ i) {
        meanR_m(i) = centroid_m[2 * i] * perParticle;
        meanP_m(i) = centroid_m[2 * i + 1] * perParticle;
    }

    meanTime_m = localStatistics[l++] * perParticle;

    meanKineticEnergy_m = localStatistics[l++] * perParticle;
    meanGamma_m = localStatistics[l++] * perParticle;
}

/* 2 * Dim centroids + Dim * ( 2 * Dim + 1 ) 2nd moments + 2 * Dim (3rd and 4th order moments)
 * --> 1st order moments: 0, ..., 2 * Dim - 1
 * --> 2nd order moments: 2 * Dim, ..., Dim * ( 2 * Dim + 1 )
 * --> 3rd order moments: Dim * ( 2 * Dim + 1 ) + 1, ..., Dim * ( 2 * Dim + 1 ) + Dim
 * (only, <x^3>, <y^3> and <z^3>)
 * --> 4th order moments: Dim * ( 2 * Dim + 1 ) + Dim + 1, ..., Dim * ( 2 * Dim + 1 ) + 2 * Dim
 *
 * For a 6x6 matrix we have each 2nd order moment (except diagonal
 * entries) twice. We only store the upper half of the matrix.
 */
template<class InputIt>
void DistributionMoments::computeStatistics(const InputIt &first, const InputIt &last)
{
    reset();

    computeMeans(first, last);

    double perParticle = 1.0 / totalNumParticles_m;
    std::vector<double> localStatistics(37);
    for (InputIt it = first; it != last; ++ it) {
        OpalParticle const& particle = *it;
        if (isParticleExcluded(particle)) {
            continue;
        }
        unsigned int l = 6;
        for (unsigned int i = 0; i < 6; ++ i) {
            localStatistics[i] += std::pow(particle[i] - centroid_m[i] * perParticle, 2);
            for (unsigned int j = 0; j <= i; ++ j, ++ l) {
                localStatistics[l] += particle[i] * particle[j];
            }
        }

        localStatistics[l++] += std::pow(particle.getTime() - meanTime_m, 2);

        for (unsigned int i = 0; i < 3; ++ i, l += 2) {
            double r2 = std::pow(particle[i], 2);
            localStatistics[l] += r2 * particle[i];
            localStatistics[l + 1] += r2 * r2;
        }

        double eKin = Util::getKineticEnergy(particle.getP(), particle.getMass());
        localStatistics[l++] += std::pow(eKin - meanKineticEnergy_m, 2);
        localStatistics[l++] += particle.getCharge();
        localStatistics[l++] += particle.getMass();
    }

    allreduce(localStatistics.data(), localStatistics.size(), std::plus<double>());

    fillMembers(localStatistics);
}

void DistributionMoments::fillMembers(std::vector<double> const& localMoments) {
    Vector_t squaredEps, fac, sumRP;
    double perParticle = 1.0 / totalNumParticles_m;

    unsigned int l = 0;
    for (; l < 6; l += 2) {
        stdR_m(l / 2) = std::sqrt(localMoments[l] * perParticle);
        stdP_m(l / 2) = std::sqrt(localMoments[l + 1] * perParticle);
    }

    for (unsigned int i = 0; i < 6; ++ i) {
        for (unsigned int j = 0; j <= i; ++ j, ++ l) {
            moments_m(i, j) = localMoments[l];
            moments_m(j, i) = moments_m(i, j);
        }
    }

    stdTime_m = localMoments[l++] * perParticle;

    for (unsigned int i = 0; i < 3; ++ i, l += 2) {
        double w1 = centroid_m[2 * i] * perParticle;
        double w2 = moments_m(2 * i , 2 * i) * perParticle;
        double w3 = localMoments[l] * perParticle;
        double w4 = localMoments[l + 1] * perParticle;
        double tmp = w2 - std::pow(w1, 2);

        halo_m(i) = (w4 + w1 * (-4 * w3 + 3 * w1 * (tmp + w2))) / tmp;
        halo_m(i) -= Options::haloShift;
    }

    stdKineticEnergy_m = std::sqrt(localMoments[l++] * perParticle);

    totalCharge_m = localMoments[l++];
    totalMass_m = localMoments[l++];

    for (unsigned int i = 0; i < 3; ++ i) {
        sumRP(i) = moments_m(2 * i, 2 * i + 1) * perParticle -  meanR_m(i) * meanP_m(i);
        stdRP_m(i) = sumRP(i) / (stdR_m(i) * stdP_m(i));
        squaredEps(i) = std::pow(stdR_m(i) * stdP_m(i), 2) - std::pow(sumRP(i), 2);
        normalizedEps_m(i) = std::sqrt(std::max(squaredEps(i), 0.0));
    }

    double betaGamma = std::sqrt(std::pow(meanGamma_m, 2) - 1.0);
    geometricEps_m = normalizedEps_m / Vector_t(betaGamma);
}

void DistributionMoments::computeMeanKineticEnergy(PartBunchBase<double, 3> const& bunch)
{
    double data[] = {0.0, 0.0};
    for (OpalParticle const& particle: bunch) {
        data[0] += Util::getKineticEnergy(particle.getP(), particle.getMass());
    }
    data[1] = bunch.getLocalNum();
    allreduce(data, 2, std::plus<double>());

    meanKineticEnergy_m = data[0] / data[1];
}

void DistributionMoments::reset()
{
    std::fill(std::begin(centroid_m), std::end(centroid_m), 0.0);
    meanR_m = 0.0;
    meanP_m = 0.0;
    stdR_m = 0.0;
    stdP_m = 0.0;
    stdRP_m = 0.0;
    normalizedEps_m = 0.0;
    geometricEps_m = 0.0;
    halo_m = 0.0;

    meanKineticEnergy_m = 0.0;
    stdKineticEnergy_m = 0.0;
    Dx_m = 0.0;
    DDx_m = 0.0;
    Dy_m = 0.0;
    DDy_m = 0.0;
    moments_m = FMatrix<double, 6, 6>(0.0);

    totalCharge_m = 0.0;
    totalMass_m = 0.0;
    totalNumParticles_m = 0;
}

bool DistributionMoments::isParticleExcluded(const OpalParticle &particle) const
{
    //FIXME After issue 287 is resolved this shouldn't be necessary anymore
    return !Options::amr
        && OpalData::getInstance()->isInOPALCyclMode()
        && particle.getId() == 0;
}