#ifndef DISTRIBUTIONMOMENTS_H
#define DISTRIBUTIONMOMENTS_H

#include "OpalParticle.h"

#include "FixedAlgebra/FMatrix.h"

#include <vector>

class DistributionMoments {
public:
    DistributionMoments();

    void compute(std::vector<OpalParticle> const&);

    Vector_t getMeanPosition() const;
    Vector_t getStandardDeviationPosition() const;
    Vector_t getMeanMomentum() const;
    Vector_t getStandardDeviationMomentum() const;
    Vector_t getNormalizedEmittance() const;
    Vector_t getGeometricEmittance() const;
    double getMeanKineticEnergy() const;

private:
    void computeMoments(std::vector<OpalParticle> const&);
    void reset();

    Vector_t meanR_m;
    Vector_t meanP_m;
    Vector_t stdR_m;
    Vector_t stdP_m;
    Vector_t normalizedEps_m;
    Vector_t geometricEps_m;
    Vector_t halo_m;

    double meanKineticEnergy_m;
    double meanGamma_m;
    double centroid_m[6];
    FMatrix<double, 6, 6> moments_m;

    unsigned int totalNumParticles_m;
};

inline
Vector_t DistributionMoments::getMeanPosition() const
{
    return meanR_m;
}

inline
Vector_t DistributionMoments::getStandardDeviationPosition() const
{
    return stdR_m;
}

inline
Vector_t DistributionMoments::getMeanMomentum() const
{
    return meanP_m;
}

inline
Vector_t DistributionMoments::getStandardDeviationMomentum() const
{
    return stdP_m;
}

inline
Vector_t DistributionMoments::getNormalizedEmittance() const
{
    return normalizedEps_m;
}

inline
Vector_t DistributionMoments::getGeometricEmittance() const
{
    return geometricEps_m;
}

inline
double DistributionMoments::getMeanKineticEnergy() const
{
    return meanKineticEnergy_m;
}

#endif