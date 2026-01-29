// ------------------------------------------------------------------------
// $RCSfile: Flagger.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Flagger
//   This class sets or resets all selection flags in a USE object.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/Flagger.h"
#include "Beamlines/FlaggedElmPtr.h"


// Class Flagger
// ------------------------------------------------------------------------

Flagger::Flagger(const Beamline &beamline, bool set):
    DefaultVisitor(beamline, false, false), flag(set)
{}


Flagger::~Flagger()
{}


void Flagger::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
    fep.setSelectionFlag(flag);
    DefaultVisitor::visitFlaggedElmPtr(fep);
}
