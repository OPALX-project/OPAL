#ifndef DISTRIBUTIONMOMENTS_H
#define DISTRIBUTIONMOMENTS_H

#include "FixedAlgebra/FMatrix.h"

#include "Vektor.h"

#include <vector>

class OpalParticle;
template<class T, unsigned Dim>
class PartBunchBase;

class DistributionMoments {
public:
    DistributionMoments();

    void compute(std::vector<OpalParticle> const&);
    void compute(PartBunchBase<double, 3> const&);
    void computeMeanKineticEnergy(PartBunchBase<double, 3> const&);

    Vector_t getMeanPosition() const;
    Vector_t getStandardDeviationPosition() const;
    Vector_t getMeanMomentum() const;
    Vector_t getStandardDeviationMomentum() const;
    Vector_t getNormalizedEmittance() const;
    Vector_t getGeometricEmittance() const;
    Vector_t getStandardDeviationRP() const;
    Vector_t getHalo() const;
    double getMeanGamma() const;
    double getMeanKineticEnergy() const;
    double getStandardDeviationKineticEnergy() const;
    double getDx() const;
    double getDDx() const;
    double getDy() const;
    double getDDy() const;


private:
    template<class Container>
    void computeMoments(Container const&);
    void reset();

    Vector_t meanR_m;
    Vector_t meanP_m;
    Vector_t stdR_m;
    Vector_t stdP_m;
    Vector_t stdRP_m;
    Vector_t normalizedEps_m;
    Vector_t geometricEps_m;
    Vector_t halo_m;

    double meanKineticEnergy_m;
    double stdKineticEnergy_m;
    double meanGamma_m;
    double Dx_m;
    double DDx_m;
    double Dy_m;
    double DDy_m;
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
Vector_t DistributionMoments::getStandardDeviationRP() const
{
    return stdRP_m;
}

inline
Vector_t DistributionMoments::getHalo() const
{
    return halo_m;
}

inline
double DistributionMoments::getMeanGamma() const
{
    return meanGamma_m;
}

inline
double DistributionMoments::getMeanKineticEnergy() const
{
    return meanKineticEnergy_m;
}

inline
double DistributionMoments::getStandardDeviationKineticEnergy() const
{
    return stdKineticEnergy_m;
}

inline
double DistributionMoments::getDx() const
{
    return Dx_m;
}

inline
double DistributionMoments::getDDx() const
{
    return DDx_m;
}

inline
double DistributionMoments::getDy() const
{
    return Dy_m;
}

inline
double DistributionMoments::getDDy() const
{
    return DDy_m;
}

#endif