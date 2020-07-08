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
#ifndef OPAL_RangeSelector_HH
#define OPAL_RangeSelector_HH

#include "Algorithms/DefaultVisitor.h"
#include "AbstractObjects/RangeRep.h"

class Beamline;


class RangeSelector: public DefaultVisitor {

public:

    /// Constructor.
    //  Attach visitor to a beamline, remember the range.
    RangeSelector(const Beamline &, const RangeRep &range);

    virtual ~RangeSelector();

    /// Execute the algorithm.
    virtual void execute();

    /// Apply the visitor to an FlaggedElmPtr.
    virtual void visitFlaggedElmPtr(const FlaggedElmPtr &);

protected:

    /// The operation to be done for beamlines.
    //  When overriding, make sure the beamline members are handled.
    virtual void handleBeamline(const FlaggedElmPtr &);

    /// The operation to be done for elements.
    virtual void handleElement(const FlaggedElmPtr &);

    /// Working data for range.
    RangeRep itsRange;

private:

    // Not implemented.
    RangeSelector();
    RangeSelector(const RangeSelector &);
    void operator=(const RangeSelector &);
};

#endif // OPAL_RangeSelector_HH
