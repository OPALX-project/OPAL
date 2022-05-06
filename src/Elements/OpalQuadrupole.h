//
// Class OpalQuadrupole
//   The QUADRUPOLE element.
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
#ifndef OPAL_OpalQuadrupole_HH
#define OPAL_OpalQuadrupole_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

class OpalQuadrupole: public OpalElement {

public:

    /// The attributes of class OpalQuadrupole.
    enum {
        K1 = COMMON,  // The normal quadrupole coefficient.
        DK1,          // The normal quadupole coefficient error.
        K1S,          // The skew quadrupole coefficient.
        DK1S,         // The skew quadrupole coefficient error.
        NSLICES,      // The number of slices / steps per element for map tracking
        SIZE
    };

    /// Exemplar constructor.
    OpalQuadrupole();

    virtual ~OpalQuadrupole();

    /// Make clone.
    virtual OpalQuadrupole* clone(const std::string& name);

    /// Print the quadrupole.
    //  Handle printing in OPAL-8 format.
    virtual void print(std::ostream&) const;

    /// Update the embedded CLASSIC multipole.
    virtual void update();

private:

    // Not implemented.
    OpalQuadrupole(const OpalQuadrupole&);
    void operator=(const OpalQuadrupole&);

    // Clone constructor.
    OpalQuadrupole(const std::string& name, OpalQuadrupole* parent);

    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalQuadrupole_HH