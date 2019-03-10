// ------------------------------------------------------------------------
// $RCSfile: Undulator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Undulator
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

#include "AbsBeamline/Undulator.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"

extern Inform *gmsg;

// Class Undulator
// ------------------------------------------------------------------------

Undulator::Undulator():
    Component(),
    nSlices_m(1)
{ }


Undulator::Undulator(const Undulator &right):
    Component(right),
    nSlices_m(right.nSlices_m)
{ }


Undulator::Undulator(const std::string &name):
    Component(name),
    nSlices_m(1)
{ }


Undulator::~Undulator()
{ }


void Undulator::accept(BeamlineVisitor &visitor) const {
    visitor.visitUndulator(*this);
}

void Undulator::initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) {
    endField = startField + getElementLength();
    RefPartBunch_m = bunch;
    startField_m = startField;
}


//set the number of slices for map tracking
void Undulator::setNSlices(const std::size_t& nSlices) { 
    nSlices_m = nSlices;
}

//get the number of slices for map tracking
std::size_t Undulator::getNSlices() const {
    return nSlices_m;
}

void Undulator::finalise() {
}

bool Undulator::bends() const {
    return false;
}

void Undulator::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = startField_m + getElementLength();
}

ElementBase::ElementType Undulator::getType() const {
    return UNDULATOR;
}
