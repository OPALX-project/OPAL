// ------------------------------------------------------------------------
// $RCSfile: UndulatorRep.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: UndulatorRep
//   Defines a concrete drift space representation.
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/UndulatorRep.h"
#include "AbsBeamline/ElementImage.h"
#include "Channels/IndirectChannel.h"


// Attribute access table.
// ------------------------------------------------------------------------

namespace {
    struct Entry {
        const char *name;
        double(UndulatorRep::*get)() const;
        void (UndulatorRep::*set)(double);
    };

    const Entry entries[] = {
        {
            "L",
            &UndulatorRep::getElementLength,
            &UndulatorRep::setElementLength
        },
        { 0, 0, 0 }
    };
}


// Class UndulatorRep
// ------------------------------------------------------------------------

UndulatorRep::UndulatorRep():
    Undulator(),
    geometry(0.0),
    k_m(0.0),
    lambda_m(0.0),
    fname_m(""),
    tranTrun_m(0.0),
    longTrun_m(0.0),
    radiationDirectory_m(""),
    truncationOrder_m(2),
    spaceCharge_m(1),
    m_m(5),
    totalTime_m(0.0),
    lFringe_m(1.0),
    is_done_m(false)
{}


UndulatorRep::UndulatorRep(const UndulatorRep &right):
    Undulator(right),
    geometry(right.geometry),
    k_m(right.k_m),
    lambda_m(right.lambda_m),
    fname_m(right.fname_m),
    tranTrun_m(right.tranTrun_m),
    longTrun_m(right.longTrun_m),
    radiationZ_m(right.radiationZ_m),
    radiationLambda_m(right.radiationLambda_m),
    radiationDirectory_m(right.radiationDirectory_m),
    meshLength_m(right.meshLength_m),
    meshResolution_m(right.meshResolution_m),
    truncationOrder_m(right.truncationOrder_m),
    spaceCharge_m(right.spaceCharge_m),
    m_m(right.m_m),
    totalTime_m(right.totalTime_m),
    lFringe_m(right.lFringe_m),
    is_done_m(right.is_done_m)
{}


UndulatorRep::UndulatorRep(const std::string &name):
    Undulator(name),
    geometry(0.0),
    k_m(0.0),
    lambda_m(0.0),
    fname_m(""),
    tranTrun_m(0.0),
    longTrun_m(0.0),
    radiationDirectory_m(""),
    truncationOrder_m(2),
    spaceCharge_m(1),
    m_m(5),
    totalTime_m(0.0),
    lFringe_m(1.0),
    is_done_m(false)
{}


UndulatorRep::~UndulatorRep()
{}


ElementBase *UndulatorRep::clone() const {
    return new UndulatorRep(*this);
}


Channel *UndulatorRep::getChannel(const std::string &aKey, bool create) {
    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        if(aKey == entry->name) {
            return new IndirectChannel<UndulatorRep>(*this, entry->get, entry->set);
        }
    }

    return ElementBase::getChannel(aKey, create);
}


NullField &UndulatorRep::getField() {
    return field;
}

const NullField &UndulatorRep::getField() const {
    return field;
}


StraightGeometry &UndulatorRep::getGeometry() {
    return geometry;
}

const StraightGeometry &UndulatorRep::getGeometry() const {
    return geometry;
}


ElementImage *UndulatorRep::getImage() const {
    ElementImage *image = ElementBase::getImage();

    for(const Entry *entry = entries; entry->name != 0; ++entry) {
        image->setAttribute(entry->name, (this->*(entry->get))());
    }

    return image;
}

void UndulatorRep::setK(double k) { k_m = k; }
double UndulatorRep::getK() const { return k_m; }

void UndulatorRep::setLambda(double lambda) { lambda_m = lambda; }
double UndulatorRep::getLambda() const { return lambda_m; }

void UndulatorRep::setFilename(const std::string&fname) { fname_m = fname; }
const std::string& UndulatorRep::getFilename() const { return fname_m;}

void UndulatorRep::setTranTrun(double tranTrun) { tranTrun_m = tranTrun; }
double UndulatorRep::getTranTrun() const { return tranTrun_m; };

void UndulatorRep::setLongTrun(double longTrun) {longTrun_m = longTrun; }
double UndulatorRep::getLongTrun() const { return longTrun_m; }

void UndulatorRep::setRadiationZ(std::vector<double> rz) { radiationZ_m = rz; }
std::vector<double> UndulatorRep::getRadiationZ() const { return radiationZ_m; }

void UndulatorRep::setRadiationLambda(std::vector<double> rl) { radiationLambda_m = rl; }
std::vector<double> UndulatorRep::getRadiationLambda() const { return radiationLambda_m; }

void UndulatorRep::setRadiationDirectory(const std::string& rd) { radiationDirectory_m = rd; }
const std::string& UndulatorRep::getRadiationDirectory() const { return radiationDirectory_m; }

void UndulatorRep::setMeshLength(std::vector<double> ml) { meshLength_m = ml; }
std::vector<double> UndulatorRep::getMeshLength() const { return meshLength_m; }

void UndulatorRep::setMeshResolution(std::vector<double> mr) { meshResolution_m = mr; }
std::vector<double> UndulatorRep::getMeshResolution() const { return meshResolution_m; }

void UndulatorRep::setTruncationOrder(unsigned int trunOrder) { truncationOrder_m = trunOrder; } 
unsigned int UndulatorRep::getTruncationOrder() const { return truncationOrder_m; }

void UndulatorRep::setSpaceCharge(bool sc) { spaceCharge_m = sc; }
bool UndulatorRep::getSpaceCharge() const { return spaceCharge_m;  }

void UndulatorRep::setTimeStepRatio(double m) { m_m = m; }
unsigned int UndulatorRep::getTimeStepRatio() const { return m_m; }

void UndulatorRep::setTotalTime(double tt) { totalTime_m = tt; }
double UndulatorRep::getTotalTime() const { return totalTime_m; }

void UndulatorRep::setLFringe(double lf) { lFringe_m = lf; }
double UndulatorRep::getLFringe() const { return lFringe_m; }

void UndulatorRep::setIsDone() { is_done_m = true; }
bool UndulatorRep::getIsDone() const { return is_done_m; }
