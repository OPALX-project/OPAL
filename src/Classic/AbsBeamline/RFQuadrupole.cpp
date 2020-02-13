// ------------------------------------------------------------------------
// $RCSfile: RFQuadrupole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RFQuadrupole
//   Defines the abstract interface for a RF quadrupole.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RFQuadrupole.h"
#include "AbsBeamline/BeamlineVisitor.h"

extern Inform *gmsg;

// Class RFQuadrupole
// ------------------------------------------------------------------------

RFQuadrupole::RFQuadrupole():
    Component()
{}


RFQuadrupole::RFQuadrupole(const RFQuadrupole &rhs):
    Component(rhs)
{}


RFQuadrupole::RFQuadrupole(const std::string &name):
    Component(name)
{}


RFQuadrupole::~RFQuadrupole()
{}


void RFQuadrupole::accept(BeamlineVisitor &visitor) const {
    visitor.visitRFQuadrupole(*this);
}

void RFQuadrupole::initialise(PartBunchBase<double, 3> *bunch, double &/*startField*/, double &/*endField*/) {
    RefPartBunch_m = bunch;
}

void RFQuadrupole::finalise()
{}

bool RFQuadrupole::bends() const {
    return false;
}


void RFQuadrupole::getDimensions(double &/*zBegin*/, double &/*zEnd*/) const {

}


ElementBase::ElementType RFQuadrupole::getType() const {
    return RFQUADRUPOLE;
}
