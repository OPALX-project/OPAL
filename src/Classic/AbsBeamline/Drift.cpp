// ------------------------------------------------------------------------
// $RCSfile: Drift.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Drift
//   Defines the abstract interface for a drift space.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Drift.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"

extern Inform *gmsg;

// Class Drift
// ------------------------------------------------------------------------

Drift::Drift():
    Drift("")
{ }


Drift::Drift(const Drift &right):
    Component(right),
    nSlices_m(right.nSlices_m)
{ }


Drift::Drift(const std::string &name):
    Component(name),
    nSlices_m(1)
{ }


Drift::~Drift()
{ }


void Drift::accept(BeamlineVisitor &visitor) const {
    visitor.visitDrift(*this);
}

void Drift::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    endField = startField + getElementLength();
    RefPartBunch_m = bunch;
    startField_m = startField;
}


//set the number of slices for map tracking
void Drift::setNSlices(const std::size_t& nSlices) { 
    nSlices_m = nSlices;
}

//get the number of slices for map tracking
std::size_t Drift::getNSlices() const {
    return nSlices_m;
}

void Drift::finalise() {
}

bool Drift::bends() const {
    return false;
}

void Drift::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + getElementLength();
}

ElementBase::ElementType Drift::getType() const {
    return DRIFT;
}
