//
// Class RangeSelector
//   This abstract class runs through a beam line and calls the pure
//   virtual methods RangeSelector::handleXXX() for each element or
//   beamline in a range.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#include "Tables/RangeSelector.h"
#include "AbsBeamline/ElementBase.h"
#include "AbstractObjects/Element.h"
#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "Beamlines/Beamline.h"
#include "Beamlines/FlaggedElmPtr.h"

class Element;


RangeSelector::RangeSelector(const Beamline &beamline, const RangeRep &range):
    DefaultVisitor(beamline, false, false), itsRange(range)
{}


RangeSelector::~RangeSelector()
{}


void RangeSelector::execute() {
    itsRange.initialize();
    DefaultVisitor::execute();
}


void RangeSelector::visitFlaggedElmPtr(const FlaggedElmPtr &fep) {
    // Are we past the beginning of the range ?
    itsRange.enter(fep);

    // Do the required operations on the beamline or element.
    ElementBase *base = fep.getElement();

    if(dynamic_cast<Beamline *>(base)) {
        handleBeamline(fep);
    } else {
        handleElement(fep);
    }

    // Are we past the end of the range ?
    itsRange.leave(fep);
}


void RangeSelector::handleBeamline(const FlaggedElmPtr &fep) {
    DefaultVisitor::visitFlaggedElmPtr(fep);
}


void RangeSelector::handleElement(const FlaggedElmPtr &fep) {
    // Default: delegate algorithm to the element, if in range.
    if(itsRange.isActive()) {
        DefaultVisitor::visitFlaggedElmPtr(fep);
    }
}
