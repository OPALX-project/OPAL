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

    void compute(const std::vector<OpalParticle>::const_iterator &,
                 const std::vector<OpalParticle>::const_iterator &);
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
    double getMeanTime() const;
    double getStdTime() const;
    double getMeanGamma() const;
    double getMeanKineticEnergy() const;
    double getStdKineticEnergy() const;
    double getDx() const;
    double getDDx() const;
    double getDy() const;
    double getDDy() const;
    double getTotalCharge() const;
    double getTotalMass() const;
    double getTotalNumParticles() const;
private:
    bool isParticleExcluded(const OpalParticle &) const;
    template<class InputIt>
    void computeMeans(const InputIt &, const InputIt &);
    template<class InputIt>
    void computeStatistics(const InputIt &, const InputIt &);
    void fillMembers(std::vector<double> const&);
    void reset();

    Vector_t meanR_m;
    Vector_t meanP_m;
    Vector_t stdR_m;
    Vector_t stdP_m;
    Vector_t stdRP_m;
    Vector_t normalizedEps_m;
    Vector_t geometricEps_m;
    Vector_t halo_m;

    double meanTime_m;
    double stdTime_m;
    double meanKineticEnergy_m;
    double stdKineticEnergy_m;
    double meanGamma_m;
    double Dx_m;
    double DDx_m;
    double Dy_m;
    double DDy_m;
    double centroid_m[6];
    FMatrix<double, 6, 6> moments_m;

    double totalCharge_m;
    double totalMass_m;
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
double DistributionMoments::getMeanTime() const
{
    return meanTime_m;
}

inline
double DistributionMoments::getStdTime() const
{
    return stdTime_m;
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
double DistributionMoments::getStdKineticEnergy() const
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

inline
double DistributionMoments::getTotalCharge() const
{
    return totalCharge_m;
}

inline
double DistributionMoments::getTotalMass() const
{
    return totalMass_m;
}

inline
double DistributionMoments::getTotalNumParticles() const
{
    return totalNumParticles_m;
}

#endif