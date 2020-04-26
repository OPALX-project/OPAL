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
    numPeriods_m(0.0),
    fname_m(""),
    truncationOrder_m(2),
    totalTime_m(0.0)
{}


UndulatorRep::UndulatorRep(const UndulatorRep &right):
    Undulator(right),
    geometry(right.geometry),
    k_m(right.k_m),
    lambda_m(right.lambda_m),
    numPeriods_m(right.numPeriods_m),
    fname_m(right.fname_m),
    meshLength_m(right.meshLength_m),
    meshResolution_m(right.meshResolution_m),
    truncationOrder_m(right.truncationOrder_m),
    totalTime_m(right.totalTime_m)
{}


UndulatorRep::UndulatorRep(const std::string &name):
    Undulator(name),
    geometry(0.0),
    k_m(0.0),
    lambda_m(0.0),
    numPeriods_m(0.0),
    fname_m(""),
    truncationOrder_m(2),
    totalTime_m(0.0)
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

void UndulatorRep::setNumPeriods(int np) { numPeriods_m = np; }
int UndulatorRep::getNumPeriods() const { return numPeriods_m; }

void UndulatorRep::setFilename(const std::string&fname) { fname_m = fname; }
const std::string& UndulatorRep::getFilename() const { return fname_m;}

void UndulatorRep::setMeshLength(std::vector<double> ml) { meshLength_m = ml; }
std::vector<double> UndulatorRep::getMeshLength() const { return meshLength_m; }

void UndulatorRep::setMeshResolution(std::vector<double> mr) { meshResolution_m = mr; }
std::vector<double> UndulatorRep::getMeshResolution() const { return meshResolution_m; }

void UndulatorRep::setTruncationOrder(unsigned int trunOrder) { truncationOrder_m = trunOrder; } 
unsigned int UndulatorRep::getTruncationOrder() const { return truncationOrder_m; }

void UndulatorRep::setTotalTime(double tt) { totalTime_m = tt; }
double UndulatorRep::getTotalTime() const { return totalTime_m; }
