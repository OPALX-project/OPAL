//
// Class OpalDrift
//   The class of OPAL drift spaces.
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
#ifndef OPAL_OpalDrift_HH
#define OPAL_OpalDrift_HH

#include "Elements/OpalElement.h"

class BoundaryGeometry;

class OpalWake;
class ParticleMatterInteraction;

class OpalDrift: public OpalElement {

public:

    enum {
         GEOMETRY = COMMON,  // geometry of boundary, one more enum member besides the common ones in OpalElement.
         NSLICES,            // The number of slices / steps per element for map tracking
         SIZE
    };

    /// Exemplar constructor.
    OpalDrift();

    virtual ~OpalDrift();

    /// Make clone.
    virtual OpalDrift* clone(const std::string& name);

    /// Test for drift.
    //  Return true.
    virtual bool isDrift() const;

    /// Update the embedded CLASSIC drift.
    virtual void update();

private:

    // Not implemented.
    OpalDrift(const OpalDrift&);
    void operator=(const OpalDrift&);

    // Clone constructor.
    OpalDrift(const std::string& name, OpalDrift* parent);

    OpalWake* owk_m;
    ParticleMatterInteraction* parmatint_m;
    BoundaryGeometry* obgeo_m;
};

#endif // OPAL_OpalDrift_HH