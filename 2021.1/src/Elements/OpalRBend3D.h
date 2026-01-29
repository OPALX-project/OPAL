//
// Class OpalRBend
//   The parent class of all OPAL bending magnets.
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
#ifndef OPAL_OpalRBend3D_HH
#define OPAL_OpalRBend3D_HH

#include "Elements/OpalBend.h"

class OpalWake;
class ParticleMatterInteraction;

class OpalRBend3D: public OpalElement {

public:

    enum {
        ANGLE = COMMON,   // The bend angle.
        K0, K0S,          // The multipole coefficients; must be in this order.
        E1,               // The edge angles.
        FMAPFN,           // File name containing on-axis field.
        GAP,              // Full gap of magnet.
        HAPERT,           // Horizontal aperture of magnet.
        DESIGNENERGY,     // the design energy of the particles
        SIZE              // Total number of attributes.
    };

    /// Exemplar constructor.
    OpalRBend3D();

    virtual ~OpalRBend3D();

    /// Make clone.
    virtual OpalRBend3D* clone(const std::string& name);

    /// Update the embedded CLASSIC bend.
    virtual void update();

    virtual void print(std::ostream&) const;

private:

    // Not implemented.
    OpalRBend3D(const OpalRBend3D&);
    void operator=(const OpalRBend3D&);

    // Clone constructor.
    OpalRBend3D(const std::string& name, OpalRBend3D* parent);

    OpalWake* owk_m;
    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalRBend3D_HH