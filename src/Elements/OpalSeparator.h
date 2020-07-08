//
// Class OpalSeparator
//   The ELSEPARATOR element.
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
#ifndef OPAL_OpalSeparator_HH
#define OPAL_OpalSeparator_HH

#include "Elements/OpalElement.h"


class OpalSeparator: public OpalElement {

public:

    /// The attributes of class OpalSeparator.
    enum {
        EX = COMMON,  // The horizontal field.
        EY,           // The vertical field.
        SIZE
    };

    /// Exemplar constructor.
    OpalSeparator();

    virtual ~OpalSeparator();

    /// Make clone.
    virtual OpalSeparator *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &);

    /// Update the embedded CLASSIC separator.
    virtual void update();

private:

    // Not implemented.
    OpalSeparator(const OpalSeparator &);
    void operator=(const OpalSeparator &);

    // Clone constructor.
    OpalSeparator(const std::string &name, OpalSeparator *parent);
};

#endif // OPAL_OpalSeparator_HH
