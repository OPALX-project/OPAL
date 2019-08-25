#ifndef OPAL_H5PARTWRAPPERFORPC_H
#define OPAL_H5PARTWRAPPERFORPC_H

/*!
  H5PartWrapperForPC: a class that manages calls to H5Part for the cyclotron tracker
*/

#include "Structure/H5PartWrapper.h"

#include "Algorithms/PBunchDefs.h"
#include "Utilities/OpalException.h"

#include "H5hut.h"

class H5PartWrapperForPC: public H5PartWrapper {
public:
    H5PartWrapperForPC(const std::string &fileName, h5_int32_t flags = H5_O_WRONLY);
    H5PartWrapperForPC(const std::string &fileName, int restartStep, std::string sourceFile, h5_int32_t flags = H5_O_RDWR);
    virtual ~H5PartWrapperForPC();

    virtual void readHeader();
    virtual void readStep(PartBunchBase<double, 3>*, h5_ssize_t firstParticle, h5_ssize_t lastParticle);

    virtual void writeHeader();
    virtual void writeStep(PartBunchBase<double, 3>*, const std::map<std::string, double> &additionalStepAttributes);

    virtual bool predecessorIsSameFlavour() const;

    bool getPreviousH5Local() const;
    double getReferenceR() const;
    double getReferenceT() const;
    double getReferenceZ() const;
    double getReferencePr() const;
    double getReferencePt() const;
    double getReferencePz() const;
    double getMeanKineticEnergy() const;
    double getMeanMomentum() const;
    double getAzimuth() const;
    double getElevation() const;
private:
    void readStepHeader(PartBunchBase<double, 3>*);
    void readStepData(PartBunchBase<double, 3>*, h5_ssize_t , h5_ssize_t);

    void writeStepHeader(PartBunchBase<double, 3>*, const std::map<std::string, double> &);
    void writeStepData(PartBunchBase<double, 3>*);

    bool previousH5Local_m;
    Vector_t referenceMomentum_m;
    Vector_t referenceLocation_m;
    h5_float64_t meanE_m;
    h5_float64_t meanMomentum_m;
    h5_float64_t azimuth_m;
    h5_float64_t elevation_m;
};

inline
bool H5PartWrapperForPC::predecessorIsSameFlavour() const {
    return (predecessorOPALFlavour_m == "opal-cycl");
}

inline
bool H5PartWrapperForPC::getPreviousH5Local() const {
    return previousH5Local_m;
}

inline
double H5PartWrapperForPC::getReferenceR() const {
    return referenceLocation_m[0];
}

inline
double H5PartWrapperForPC::getReferenceT() const {
    return referenceLocation_m[1];
}

inline
double H5PartWrapperForPC::getReferenceZ() const {
    return referenceLocation_m[2];
}

inline
double H5PartWrapperForPC::getReferencePr() const {
    return referenceMomentum_m[0];
}

inline
double H5PartWrapperForPC::getReferencePt() const {
    return referenceMomentum_m[1];
}

inline
double H5PartWrapperForPC::getReferencePz() const {
    return referenceMomentum_m[2];
}

inline
double H5PartWrapperForPC::getMeanKineticEnergy() const {
    return meanE_m;
}

inline
double H5PartWrapperForPC::getMeanMomentum() const {
    return meanMomentum_m;
}

inline
double H5PartWrapperForPC::getAzimuth() const {
    return azimuth_m;
}

inline
double H5PartWrapperForPC::getElevation() const {
    return elevation_m;
}

#endif //OPAL_H5PARTWRAPPERFORPC_H