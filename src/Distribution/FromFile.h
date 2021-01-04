#ifndef FROMFILE_H
#define FROMFILE_H

#include "DistributionBase.h"

class FromFile: public DistributionBase
{
 public:
    FromFile();
    std::vector<ParticlePhaseSpace> create(size_t numberOfParticles,
                                           double massInEv,
                                           double charge) override;

    Inform& print(Inform& os) const override;
    DistributionType getType() const override;
    double getEmissionTime() const override;
};

inline
DistributionBase::DistributionType FromFile::getType() const
{
    return FROMFILE;
}

#endif