// ------------------------------------------------------------------------
// $RCSfile: Patch.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Patch
//   Defines the abstract interface for a geometry patch.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Patch.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Patch
// ------------------------------------------------------------------------

Patch::Patch():
    Component()
{}


Patch::Patch(const Patch &rhs):
    Component(rhs)
{}


Patch::Patch(const std::string &name):
    Component(name)
{}


Patch::~Patch()
{}


void Patch::accept(BeamlineVisitor &visitor) const {
    visitor.visitPatch(*this);
}

void Patch::initialise(PartBunchBase<double, 3> *bunch, double &/*startField*/, double &/*endField*/) {
    RefPartBunch_m = bunch;
}

void Patch::finalise()
{}

bool Patch::bends() const {
    return false;
}

void Patch::getDimensions(double &/*zBegin*/, double &/*zEnd*/) const {

}


ElementBase::ElementType Patch::getType() const {
    return PATCH;
}
