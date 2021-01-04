#ifndef DISTRIBUTIONBASE_H
#define DISTRIBUTIONBASE_H

#include "Algorithms/PBunchDefs.h"

#include <vector>

class DistributionBase
{
 public:
    DistributionBase() {}
{
    enum class DistributionType {
        NODIST,
        FROMFILE,
        GAUSS,
        BINOMIAL,
        FLATTOP,
        MULTIGAUSS,
        GUNGAUSSFLATTOPTH,
        ASTRAFLATTOPTH,
        MATCHEDGAUSS
    };

    struct ParticlePhaseSpace
    {
        Vector_t R_m;
        Vector_t P_m;
    };
    
    virtual std::vector<ParticlePhaseSpace> create(size_t numberOfParticles,
                                                   double massInEv,
                                                   double charge) = 0;
    virtual Inform& void print(Inform& os) const = 0;
    virtual DistributionType getType() const = 0;
    virtual double getEmissionTime() const = 0;
    
};

#endif