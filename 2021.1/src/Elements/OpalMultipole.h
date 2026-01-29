//
// Class OpalMultipole
//   The MULTIPOLE element.
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
#ifndef OPAL_OpalMultipole_HH
#define OPAL_OpalMultipole_HH

#include "Elements/OpalElement.h"


class OpalMultipole: public OpalElement {

public:

    /// The attributes of class OpalMultipole.
    enum {
        KN = COMMON,  // The normal field components.
        DKN,          // The normal field component errors.
        KS,           // The skewed field components.
        DKS,          // The skewed field component errors.
        SIZE
    };

    /// Exemplar constructor.
    OpalMultipole();

    virtual ~OpalMultipole();

    /// Make clone.
    virtual OpalMultipole *clone(const std::string &name);

    /// Print the object.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream &) const;

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalMultipole(const OpalMultipole &);
    void operator=(const OpalMultipole &);

    // Clone constructor.
    OpalMultipole(const std::string &name, OpalMultipole *parent);
};

#endif // OPAL_OpalMultipole_HH
