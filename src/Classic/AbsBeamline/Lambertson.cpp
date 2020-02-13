// ------------------------------------------------------------------------
// $RCSfile: Lambertson.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Lambertson
//   Defines the abstract interface for a Lambertson septum magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Lambertson.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Lambertson
// ------------------------------------------------------------------------

Lambertson::Lambertson():
    Component()
{ }


Lambertson::Lambertson(const Lambertson &rhs):
    Component(rhs)
{ }


Lambertson::Lambertson(const std::string &name):
    Component(name)
{ }


Lambertson::~Lambertson()
{ }


void Lambertson::accept(BeamlineVisitor &visitor) const {
    visitor.visitLambertson(*this);
}

void Lambertson::initialise(PartBunchBase<double, 3> *bunch, double &/*startField*/, double &/*endField*/) {
    RefPartBunch_m = bunch;
}

void Lambertson::finalise()
{ }

bool Lambertson::bends() const {
    return false;
}

void Lambertson::getDimensions(double &/*zBegin*/, double &/*zEnd*/) const {
}

ElementBase::ElementType Lambertson::getType() const {
    return LAMBERTSON;
}
