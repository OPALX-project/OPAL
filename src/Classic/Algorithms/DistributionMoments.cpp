#include "DistributionMoments.h"

#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include "Message/GlobalComm.h"

DistributionMoments::DistributionMoments():
    meanR_m(0.0),
    meanP_m(0.0),
    stdR_m(0.0),
    stdP_m(0.0),
    normalizedEps_m(0.0),
    geometricEps_m(0.0),
    meanKineticEnergy_m(0.0),
    moments_m(),
    totalNumParticles_m(0)
{
    std::fill(std::begin(centroid_m), std::end(centroid_m), 0.0);
}

void DistributionMoments::compute(std::vector<OpalParticle> const& particles)
{
    reset();

    unsigned int localNumParticles = particles.size();
    allreduce(&localNumParticles, &totalNumParticles_m, 1, std::plus<unsigned int>());

    if (totalNumParticles_m == 0) {
        return;
    }

    computeMoments(particles);

    Vector_t squaredEps, fac, squaredSumR, squaredSumP, sumRP;
    double perParticle = 1.0 / totalNumParticles_m;
    for (unsigned int i = 0; i < 3; ++ i) {
        meanR_m(i) = centroid_m[2 * i] * perParticle;
        meanP_m(i) = centroid_m[2 * i + 1] * perParticle;
        squaredSumR(i) = moments_m(2 * i, 2 * i) - totalNumParticles_m * std::pow(meanR_m(i), 2);
        squaredSumP(i) = std::max(0.0,
                                  moments_m(2 * i + 1, 2 * i + 1) - totalNumParticles_m * std::pow(meanP_m(i), 2));
        sumRP(i) = moments_m(2 * i, 2 * i + 1) - totalNumParticles_m * meanR_m(i) * meanP_m(i);
    }

    squaredEps = (squaredSumR * squaredSumP - sumRP * sumRP) * std::pow(perParticle, 2);
    // sumRP *= perParticle;

    for (unsigned int i = 0; i < 3; ++ i) {
        stdR_m(i) = std::sqrt(squaredSumR(i) * perParticle);
        stdP_m(i) = std::sqrt(squaredSumP(i) * perParticle);
        normalizedEps_m(i) = std::sqrt(std::max(squaredEps(i), 0.0));
        double tmp = stdR_m(i) * stdP_m(i);
        fac(i) = (std::abs(tmp) < 1e-10) ? 0.0: 1.0 / tmp;
    }

    // stdRP_m = sumRP * fac;

    double betaGamma = std::sqrt(std::pow(meanGamma_m, 2) - 1.0);
    geometricEps_m = normalizedEps_m / Vector_t(betaGamma);
}

void DistributionMoments::computeMoments(std::vector<OpalParticle> const& particles)
{
    std::vector<double> localMoments(35);

    for (OpalParticle const& particle: particles) {
        unsigned int l = 6;
        for (unsigned int i = 0; i < 6; ++ i) {
            localMoments[i] += particle[i];
            for (unsigned int j = 0; j <= i; ++ j, ++ l) {
                localMoments[l] += particle[i] * particle[j];
            }
        }

        for (unsigned int i = 0; i < 3; ++ i, ++ l) {
            double r2 = std::pow(particle[i], 2);
            localMoments[l] += r2 * particle[i];
            localMoments[l + 3] += r2 * r2;
        }
        double gamma = Util::getGamma(particle.P());
        localMoments[33] = (gamma - 1.0) * particle.mass();
        localMoments[34] = gamma;
    }


    allreduce(localMoments.data(), localMoments.size(), std::plus<double>());

    for (unsigned int i = 0; i < 6; ++ i) {
        centroid_m[i] = localMoments[i];
    }

    unsigned int l = 6;
    for (unsigned int i = 0; i < 6; ++ i) {
        for (unsigned int j = 0; j <= i; ++ j, ++ l) {
            moments_m(i, j) = localMoments[l];
            moments_m(j, i) = moments_m(i, j);
        }
    }

    double perParticle = 1.0 / totalNumParticles_m;

    for (unsigned int i = 0; i < 3; ++ i) {
        double w1 = centroid_m[2 * i] * perParticle;
        double w2 = moments_m(2 * i , 2 * i) * perParticle;
        double w3 = localMoments[i + l] * perParticle;
        double w4 = localMoments[i + l + 3] * perParticle;
        double tmp = w2 - std::pow(w1, 2);

        halo_m(i) = (w4 + w1 * (-4 * w3 + 3 * w1 * (tmp + w2))) / tmp;
        halo_m(i) -= Options::haloShift;
    }

    meanKineticEnergy_m = localMoments[33] * perParticle;
    meanGamma_m = localMoments[34] * perParticle;
}

void DistributionMoments::reset()
{
    meanR_m = 0.0;
    meanP_m = 0.0;
    stdR_m = 0.0;
    stdP_m = 0.0;
    normalizedEps_m = 0.0;
    geometricEps_m = 0.0;
    halo_m = 0.0;

    meanKineticEnergy_m = 0.0;
    std::fill(std::begin(centroid_m), std::end(centroid_m), 0.0);
    moments_m = FMatrix<double, 6, 6>(0.0);

    totalNumParticles_m = 0;
}