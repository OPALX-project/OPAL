//
// Class OpalSextupole
//   The SEXTUPOLE element.
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
#ifndef OPAL_OpalSextupole_HH
#define OPAL_OpalSextupole_HH

#include "Elements/OpalElement.h"


class OpalSextupole: public OpalElement {

public:

    /// The attributes of class OpalSextupole.
    enum {
        K2 = COMMON,  // The normal sextupole strength.
        DK2,          // The normal sextupole strength error.
        K2S,          // The skew sextupole strength.
        DK2S,         // The skew sextupole strength error.
        SIZE
    };

    /// Exemplar constructor.
    OpalSextupole();

    virtual ~OpalSextupole();

    /// Make clone.
    virtual OpalSextupole *clone(const std::string &name);

    /// Print the sextupole.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalSextupole(const OpalSextupole &);
    void operator=(const OpalSextupole &);

    // Clone constructor.
    OpalSextupole(const std::string &name, OpalSextupole *parent);
};

#endif // OPAL_OpalSextupole_HH
