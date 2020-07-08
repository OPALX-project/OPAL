//
// Class Selector
//   Set selection flags for a given range in a beam line.
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
#ifndef OPAL_Selector_HH
#define OPAL_Selector_HH

#include "Tables/RangeSelector.h"
#include <string>

class Element;
class RangeRep;
class RegularExpression;

class Selector: public RangeSelector {

public:

    /// Constructor.
    //  Attach visitor to [b]bl[/b].  Remember range [b]range[/b],
    //  class name [b]cName[/b], type name [b]tName[/b], and pattern
    //  string [b]pString[/b].
    Selector(const Beamline &, const RangeRep &range, const std::string &cName,
             const std::string &tName, const std::string &pString);

    virtual ~Selector();

    /// Execute the selection.
    virtual void execute();

    /// Return the count of selected elements.
    int getCount() const;

protected:

    /// The operation to be done for elements.
    virtual void handleElement(const FlaggedElmPtr &);

private:

    // Not implemented.
    Selector();
    Selector(const Selector &);
    void operator=(const Selector &);

    // Class of "SELECT", or zero.
    const Element *itsClass;

    // Type name of "SELECT".
    const std::string itsType;

    // Pattern of "SELECT", or zero.
    const RegularExpression *itsPattern;

    // The count of selected elements.
    int itsCount;
};

#endif // OPAL_Selector_HH
