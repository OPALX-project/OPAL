
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

#include <gsl/gsl_histogram.h>

#include <boost/numeric/conversion/cast.hpp>

extern Inform* gmsg;

const double DistributionMoments::percentileOneSigmaNormalDist_m = std::erf(1 / sqrt(2));
const double DistributionMoments::percentileTwoSigmasNormalDist_m = std::erf(2 / sqrt(2));
const double DistributionMoments::percentileThreeSigmasNormalDist_m = std::erf(3 / sqrt(2));
const double DistributionMoments::percentileFourSigmasNormalDist_m = std::erf(4 / sqrt(2));

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
    std::vector<double> maxima(6, std::numeric_limits<double>::lowest());
    for (InputIt it = first; it != last; ++ it) {
        OpalParticle const& particle = *it;
        if (isParticleExcluded(particle)) {
            -- localNum;
            continue;
        }
        for (unsigned int i = 0; i < 3; ++ i) {
            maxima[2 * i] = std::max(maxima[2 * i], particle[2 * i]);
            maxima[2 * i + 1] = std::max(maxima[2 * i + 1], -particle[2 * i]); // calculates the minimum
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
    allreduce(maxima.data(), 6, std::greater<double>());
    totalNumParticles_m = localStatistics.back();

    double perParticle = 1.0 / totalNumParticles_m;
    unsigned int l = 0;

    for (; l < 6; ++ l) {
        centroid_m[l] = localStatistics[l];
    }
    for (unsigned int i = 0; i < 3; ++ i) {
        meanR_m(i) = centroid_m[2 * i] * perParticle;
        meanP_m(i) = centroid_m[2 * i + 1] * perParticle;
        maxR_m(i) = maxima[2 * i];
        minR_m(i) = -maxima[2 * i + 1];
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

    computePercentiles(first, last);
}


template<class InputIt>
void DistributionMoments::computePercentiles(const InputIt & first, const InputIt & last) {
    if (!Options::computePercentiles || totalNumParticles_m < 100) {
        return;
    }

    std::vector<gsl_histogram*> histograms(3);
    // For a normal distribution the number of exchanged data between the cores is minimized
    // if the number of histogram bins follows the following formula. Since we can't know
    // how many particles are in each bin for the real distribution we use this formula.
    unsigned int numBins = 3.5 * std::pow(3, std::log10(totalNumParticles_m));

    Vector_t maxR;
    for (unsigned int d = 0; d < 3; ++ d) {
        maxR(d) = 1.0000001 * std::max(maxR_m[d] - meanR_m[d], meanR_m[d] - minR_m[d]);
        histograms[d] = gsl_histogram_alloc(numBins);
        gsl_histogram_set_ranges_uniform(histograms[d], 0.0, maxR(d));
    }
    for (InputIt it = first; it != last; ++ it) {
        OpalParticle const& particle = *it;
        for (unsigned int d = 0; d < 3; ++ d) {
            gsl_histogram_increment(histograms[d], std::abs(particle[2 * d] - meanR_m[d]));
        }
    }

    std::vector<int> localHistogramValues(3 * (numBins + 1)), globalHistogramValues(3 * (numBins + 1));
    for (unsigned int d = 0; d < 3; ++ d) {
        int j = 0;
        size_t accumulated = 0;
        std::generate(localHistogramValues.begin() + d * (numBins + 1) + 1,
                      localHistogramValues.begin() + (d + 1) * (numBins + 1),
                      [&histograms,&d,&j,&accumulated](){
                          accumulated += gsl_histogram_get(histograms[d], j ++);
                          return accumulated;});

        gsl_histogram_free(histograms[d]);
    }

    allreduce(localHistogramValues.data(), globalHistogramValues.data(), 3 * (numBins + 1), std::plus<int>());

    int numParticles68 = boost::numeric_cast<int>(std::floor(totalNumParticles_m * percentileOneSigmaNormalDist_m + 0.5));
    int numParticles95 = boost::numeric_cast<int>(std::floor(totalNumParticles_m * percentileTwoSigmasNormalDist_m + 0.5));
    int numParticles99 = boost::numeric_cast<int>(std::floor(totalNumParticles_m * percentileThreeSigmasNormalDist_m + 0.5));
    int numParticles99_99 = boost::numeric_cast<int>(std::floor(totalNumParticles_m * percentileFourSigmasNormalDist_m + 0.5));

    for (int d = 0; d < 3; ++ d) {
        unsigned int localNum = last - first, current = 0;
        std::vector<Vektor<double, 2>> oneDPhaseSpace(localNum);
        for (InputIt it = first; it != last; ++ it, ++ current) {
            OpalParticle const& particle = *it;
            oneDPhaseSpace[current](0) = particle[2 * d];
            oneDPhaseSpace[current](1) = particle[2 * d + 1];
        }
        std::sort(oneDPhaseSpace.begin(), oneDPhaseSpace.end(),
                  [d, this](Vektor<double, 2>& left, Vektor<double, 2>& right) {
                      return std::abs(left[0] - meanR_m[d]) < std::abs(right[0] - meanR_m[d]);
                  });

        iterator_t endSixtyEight, endNinetyFive, endNinetyNine, endNinetyNine_NinetyNine;
        endSixtyEight = endNinetyFive = endNinetyNine = endNinetyNine_NinetyNine = oneDPhaseSpace.end();

        std::tie(sixtyEightPercentile_m[d], endSixtyEight) =
            determinePercentilesDetail(oneDPhaseSpace.begin(),
                                       oneDPhaseSpace.end(),
                                       globalHistogramValues,
                                       localHistogramValues,
                                       d, numParticles68);

        std::tie(ninetyFivePercentile_m[d], endNinetyFive) =
            determinePercentilesDetail(oneDPhaseSpace.begin(),
                                       oneDPhaseSpace.end(),
                                       globalHistogramValues,
                                       localHistogramValues,
                                       d, numParticles95);

        std::tie(ninetyNinePercentile_m[d], endNinetyNine) =
            determinePercentilesDetail(oneDPhaseSpace.begin(),
                                       oneDPhaseSpace.end(),
                                       globalHistogramValues,
                                       localHistogramValues,
                                       d, numParticles99);

        std::tie(ninetyNine_NinetyNinePercentile_m[d], endNinetyNine_NinetyNine) =
            determinePercentilesDetail(oneDPhaseSpace.begin(),
                                       oneDPhaseSpace.end(),
                                       globalHistogramValues,
                                       localHistogramValues,
                                       d, numParticles99_99);

        normalizedEps68Percentile_m[d] = computeNormalizedEmittance(oneDPhaseSpace.begin(), endSixtyEight);
        normalizedEps95Percentile_m[d] = computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyFive);
        normalizedEps99Percentile_m[d] = computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyNine);
        normalizedEps99_99Percentile_m[d] = computeNormalizedEmittance(oneDPhaseSpace.begin(), endNinetyNine_NinetyNine);
    }
}

/** Computes the percentile and the range of all local particles that are contained therein.
 *  In a first step the container globalAccumulatedHistogram is looped through until accumulated
 *  histogram value exceeds the required number of particles. The percentile then is between the
 *  boundaries of the last histogram bin before the loop stopped. Then all particle coordinates
 *  that are between the boundaries of this bin are communicated acros all nodes and sorted.
 *  The exact percentile is then determined by counting the n smallest coordinates such that
 *  the total number of partiles results is equal to 'numRequiredParticles'. In accordance with
 *  matlab (?) the percentile is the midpoint between the last particle within the percentile and
 *  tje first particle outside. Finally each node determines which of its particles are contained
 *  in the percentile.
 *
 *  To determine the histogram, the coordinates should not be used directly. Instead, the
 *  absolute value of the difference between a coordinate and the mean, |x - <x>|, should be used
 *  so that the percentile values are similar to the standard deviation.

 * @param begin: begin of a container containing the one dimensional phase space of all local
 *               particles.
 * @param end:   end of the container.
 * @param globalAccumulatedHistogram: container with partial sum of histogram values of position
 *                                    coordinates summed up across all nodes. The first value
 *                                    should be 0.
 * @param localAccumulatedHistogram: container with partial sum of histogram values of position
 *                                   coordinates of all local particles. The first value should be 0.
 * @param dimension: dimension of the one dimensional phase space.
 * @param numRequiredParticles: number of particles that are contained in the requested percentile.
 *                              Is determined by the total number of particles and the percentile.
 * @return: pair of percentile and iterator pointing to the element after the range contained in the
 *          percentile.
 */
std::pair<double, DistributionMoments::iterator_t>
DistributionMoments::determinePercentilesDetail(const DistributionMoments::iterator_t& begin,
                                                const DistributionMoments::iterator_t& end,
                                                const std::vector<int>& globalAccumulatedHistogram,
                                                const std::vector<int>& localAccumulatedHistogram,
                                                unsigned int dimension,
                                                int numRequiredParticles) const {
    unsigned int numBins = globalAccumulatedHistogram.size() / 3;
    double percentile = 0.0;
    iterator_t endPercentile = end;
    for (unsigned int i = 1; i < numBins; ++ i) {
        unsigned int idx = dimension * numBins + i;
        if (globalAccumulatedHistogram[idx] > numRequiredParticles) {
            iterator_t beginBin = begin + localAccumulatedHistogram[idx - 1];
            iterator_t endBin = begin + localAccumulatedHistogram[idx];
            unsigned int numMissingParticles = numRequiredParticles - globalAccumulatedHistogram[idx - 1];
            unsigned int shift = 2;
            while (numMissingParticles == 0) {
                beginBin = begin + localAccumulatedHistogram[idx - shift];
                numMissingParticles = numRequiredParticles - globalAccumulatedHistogram[idx - shift];
                ++ shift;
            }

            std::vector<unsigned int> numParticlesInBin(Ippl::getNodes() + 1);
            numParticlesInBin[Ippl::myNode() + 1] = endBin - beginBin;
            allreduce(&(numParticlesInBin[1]), Ippl::getNodes(), std::plus<unsigned int>());
            std::partial_sum(numParticlesInBin.begin(), numParticlesInBin.end(), numParticlesInBin.begin());

            std::vector<double> positions(numParticlesInBin.back());
            std::transform(beginBin, endBin, positions.begin() + numParticlesInBin[Ippl::myNode()],
                           [&dimension, this](Vektor<double, 2> const& particle)
                           { return std::abs(particle[0] - meanR_m[dimension]); });
            allreduce(&(positions[0]), positions.size(), std::plus<double>());
            std::sort(positions.begin(), positions.end());

            percentile = (*(positions.begin() + numMissingParticles - 1)
                          + *(positions.begin() + numMissingParticles)) / 2;
            for (iterator_t it = beginBin; it != endBin; ++ it) {
                if (std::abs((*it)[0] - meanR_m[dimension]) > percentile) {
                    return std::make_pair(percentile, it);
                }
            }
            endPercentile = endBin;
            break;
        }
    }
    return std::make_pair(percentile, endPercentile);
}

double DistributionMoments::computeNormalizedEmittance(const DistributionMoments::iterator_t& begin,
                                                       const DistributionMoments::iterator_t& end) const {
    double localStatistics[] = {0.0, 0.0, 0.0, 0.0};
    localStatistics[0] = end - begin;
    for (iterator_t it = begin; it < end; ++ it) {
        const Vektor<double, 2>& rp = *it;
        localStatistics[1] += rp(0);
        localStatistics[2] += rp(1);
        localStatistics[3] += rp(0) * rp(1);
    }
    allreduce(&(localStatistics[0]), 4, std::plus<double>());

    double numParticles = localStatistics[0];
    double perParticle = 1 / localStatistics[0];
    double meanR = localStatistics[1] * perParticle;
    double meanP = localStatistics[2] * perParticle;
    double RP = localStatistics[3] * perParticle;

    localStatistics[0] = 0.0;
    localStatistics[1] = 0.0;
    for (iterator_t it = begin; it < end; ++ it) {
        const Vektor<double, 2>& rp = *it;
        localStatistics[0] += std::pow(rp(0) - meanR, 2);
        localStatistics[1] += std::pow(rp(1) - meanP, 2);
    }
    allreduce(&(localStatistics[0]), 2, std::plus<double>());

    double stdR = std::sqrt(localStatistics[0] / numParticles);
    double stdP = std::sqrt(localStatistics[1] / numParticles);
    double sumRP = RP -  meanR * meanP;
    double squaredEps = std::pow(stdR * stdP, 2) - std::pow(sumRP, 2);
    double normalizedEps = std::sqrt(std::max(squaredEps, 0.0));

    return normalizedEps;
}

void DistributionMoments::fillMembers(std::vector<double> const& localMoments) {
    Vector_t squaredEps, fac, sumRP;
    double perParticle = 1.0 / totalNumParticles_m;

    unsigned int l = 0;
    for (; l < 6; l += 2) {
        stdR_m(l / 2) = std::sqrt(localMoments[l] / totalNumParticles_m);
        stdP_m(l / 2) = std::sqrt(localMoments[l + 1] / totalNumParticles_m);
    }

    for (unsigned int i = 0; i < 6; ++ i) {
        for (unsigned int j = 0; j <= i; ++ j, ++ l) {
            moments_m(i, j) = localMoments[l] * perParticle;
            moments_m(j, i) = moments_m(i, j);
        }
    }

    stdTime_m = localMoments[l++] * perParticle;

    for (unsigned int i = 0; i < 3; ++ i, l += 2) {
        double w1 = centroid_m[2 * i] * perParticle;
        double w2 = moments_m(2 * i , 2 * i);
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
        sumRP(i) = moments_m(2 * i, 2 * i + 1) -  meanR_m(i) * meanP_m(i);
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
    moments_m = FMatrix<double, 6, 6>(0.0);

    totalCharge_m = 0.0;
    totalMass_m = 0.0;
    totalNumParticles_m = 0;

    sixtyEightPercentile_m = 0.0;
    normalizedEps68Percentile_m = 0.0;
    ninetyFivePercentile_m = 0.0;
    normalizedEps95Percentile_m = 0.0;
    ninetyNinePercentile_m = 0.0;
    normalizedEps99Percentile_m = 0.0;
    ninetyNine_NinetyNinePercentile_m = 0.0;
    normalizedEps99_99Percentile_m = 0.0;
}

bool DistributionMoments::isParticleExcluded(const OpalParticle &particle) const
{
    //FIXME After issue 287 is resolved this shouldn't be necessary anymore
    return !Options::amr
        && OpalData::getInstance()->isInOPALCyclMode()
        && particle.getId() == 0;
}