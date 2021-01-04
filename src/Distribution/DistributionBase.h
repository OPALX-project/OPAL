#ifndef DISTRIBUTIONBASE_H
#define DISTRIBUTIONBASE_H

#include "Algorithms/PBunchDefs.h"

class DistributionBase
{
 public:
    DistributionBase() {}

    struct ParticlePhaseSpace
    {
        Vector_t R_m;
        Vector_t P_m;
    };
    
    virtual std::vector<ParticlePhaseSpace> create(size_t numberOfParticles,
                                                   double massInEv,
                                                   double charge) = 0;
    virtual void print(Inform &os) const = 0;
    virtual double getEmissionTime() const = 0;
    
};

#endif