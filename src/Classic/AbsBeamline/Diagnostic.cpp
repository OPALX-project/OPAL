// ------------------------------------------------------------------------
// $RCSfile: Diagnostic.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Diagnostic
//   Defines the abstract interface for a beam diagnostic.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Diagnostic.h"
#include "AbsBeamline/BeamlineVisitor.h"


// Class Diagnostic
// ------------------------------------------------------------------------

Diagnostic::Diagnostic():
    Component()
{ }


Diagnostic::Diagnostic(const Diagnostic &rhs):
    Component(rhs)
{ }


Diagnostic::Diagnostic(const std::string &name):
    Component(name)
{ }


Diagnostic::~Diagnostic()
{ }


void Diagnostic::accept(BeamlineVisitor &visitor) const {
    visitor.visitDiagnostic(*this);
}

void Diagnostic::initialise(PartBunchBase<double, 3> *bunch, double &/*startField*/, double &/*endField*/) {
    RefPartBunch_m = bunch;
}

void Diagnostic::finalise()
{ }

bool Diagnostic::bends() const {
    return false;
}

ElementBase::ElementType Diagnostic::getType() const {
    return DIAGNOSTIC;
}

void Diagnostic::getDimensions(double &/*zBegin*/, double &/*zEnd*/) const
{ }
