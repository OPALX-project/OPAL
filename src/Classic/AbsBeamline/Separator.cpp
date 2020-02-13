// ------------------------------------------------------------------------
// $RCSfile: Separator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Separator
//   Defines the abstract interface for an  separator.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Separator.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Separator
// ------------------------------------------------------------------------

Separator::Separator():
    Component()
{}


Separator::Separator(const Separator &right):
    Component(right)
{}


Separator::Separator(const std::string &name):
    Component(name)
{}


Separator::~Separator()
{}


void Separator::accept(BeamlineVisitor &visitor) const {
    visitor.visitSeparator(*this);
}

void Separator::initialise(PartBunchBase<double, 3> *bunch, double &/*startField*/, double &/*endField*/) {
    RefPartBunch_m = bunch;
}

void Separator::finalise()
{}

bool Separator::bends() const {
    return false;
}

void Separator::getDimensions(double &/*zBegin*/, double &/*zEnd*/) const {

}


ElementBase::ElementType Separator::getType() const {
    return SEPARATOR;
}
