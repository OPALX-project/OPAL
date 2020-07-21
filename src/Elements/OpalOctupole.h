//
// Class OpalOctupole
//   The OCTUPOLE element.
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
#ifndef OPAL_OpalOctupole_HH
#define OPAL_OpalOctupole_HH


#include "Elements/OpalElement.h"


class OpalOctupole: public OpalElement {

public:

    /// The attributes of class OpalOctupole.
    enum {
        K3 = COMMON,  // The normal octupole coefficient.
        DK3,          // The normal octupole coefficient error.
        K3S,          // The skew octupole coefficient.
        DK3S,          // The skew octupole coefficient error.
        SIZE
    };

    /// Exemplar constructor.
    OpalOctupole();

    virtual ~OpalOctupole();

    /// Make clone.
    virtual OpalOctupole *clone(const std::string &name);

    /// Print the element.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalOctupole(const OpalOctupole &);
    void operator=(const OpalOctupole &);

    // Clone constructor.
    OpalOctupole(const std::string &name, OpalOctupole *parent);
};

#endif // OPAL_OpalOctupole_HH
