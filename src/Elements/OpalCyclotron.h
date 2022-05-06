//
// Class OpalCyclotron
//   The OpalCyclotron element.
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
#ifndef OPAL_OpalCyclotron_HH
#define OPAL_OpalCyclotron_HH

#include "Elements/OpalElement.h"

class BoundaryGeometry;

class OpalCyclotron: public OpalElement {

public:

    /// The attributes of class OpalCyclotron.
    enum {
        GEOMETRY = COMMON,  // geometry of boundary
        CYHARMON,           // The harmonic number of the cyclotron
        SYMMETRY,           // The symetry of the field
        RINIT,              // The initial radius [mm]
        PRINIT,             // The initial radial momentum [pr/p0] []
        PHIINIT,            // The initial phase [deg]
        ZINIT,              // The initial z coordinate [mm]
        PZINIT,             // The initial vertical momentum [pz/p0] []
        RFFREQ,             // First hamonic of the RF system [MHz]
        FMAPFN,             // The filename of the mid-plane fieldmap
        RFMAPFN,            // The filename(s) of the RF fieldmap
        RFFCFN,             // The filename(s) of coefficients for RF frequency function f(t)
        RFVCFN,             // The filename(s) of coefficients for RF voltage function v(t)
        BSCALE,             // A scalar to scale the B-field
        ESCALE,             // A scalar to scale the RF field
        RFPHI,              // the initial phase of RF field
        SUPERPOSE,          // whether the electric field map are superposed or not
        MINZ,               // minimal vertical extend of the machine
        MAXZ,               // maximal vertical extend of the machine
        MINR,               // minimal radial extend of the machine
        MAXR,               // maximal radial extend of the machine
        FMLOWE,             // minimal energy of the field map
        FMHIGHE,            // maximal energy of the field map
        SPIRAL,             // flag whether or not this is a spiral inflector simulation
        TRIMCOILTHRESHOLD,  // minimum B-field for which trim coils are applied
        TRIMCOIL,           // list of trim coils
        SIZE
    };

    OpalCyclotron();

    virtual ~OpalCyclotron();

    /// Make clone.
    virtual OpalCyclotron* clone(const std::string& name);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalCyclotron(const OpalCyclotron&);
    void operator=(const OpalCyclotron&);

    // Clone constructor.
    OpalCyclotron(const std::string& name, OpalCyclotron* parent);

    BoundaryGeometry* obgeo_m;

};

#endif // OPAL_OpalCyclotron_HH
