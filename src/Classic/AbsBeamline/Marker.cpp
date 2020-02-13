// ------------------------------------------------------------------------
// $RCSfile: Marker.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Marker
//   Defines the abstract interface for a marker element.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Marker.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Marker
// ------------------------------------------------------------------------


Marker::Marker():
    Component()
{ }


Marker::Marker(const Marker &right):
    Component(right)
{ }


Marker::Marker(const std::string &name):
    Component(name)
{ }


Marker::~Marker()
{ }


void Marker::accept(BeamlineVisitor &visitor) const {
    visitor.visitMarker(*this);
}

void Marker::initialise(PartBunchBase<double, 3> *bunch, double &/*startField*/, double &/*endField*/) {
    RefPartBunch_m = bunch;
}

void Marker::finalise()
{ }

bool Marker::bends() const {
    return false;
}


void Marker::getDimensions(double &/*zBegin*/, double &/*zEnd*/) const {
}

ElementBase::ElementType Marker::getType() const {
    return MARKER;
}

