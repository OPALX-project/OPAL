//
// Class OpalCavity
//   The RFCAVITY element.
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
#ifndef OPAL_OpalCavity_HH
#define OPAL_OpalCavity_HH

#include "Elements/OpalElement.h"

class OpalWake;
class BoundaryGeometry;

class OpalCavity: public OpalElement {

public:

    /// The attributes of class OpalCavity.
    enum {
        VOLT = COMMON,  // The peak voltage.
        DVOLT,          // The peak voltage error.
        GEOMETRY,       // geometry of boundary
        FREQ,           // The RF frequency.
        LAG,            // The phase lag.
        DLAG,           // The phase lag error.
        FMAPFN,         // The filename of the fieldmap
        FAST,           // Faster but less accurate
        APVETO,         // Do not use this cavity in the Autophase procedure
        RMIN,           // Minimal Radius
        RMAX,           // Maximal Radius
        ANGLE,          // the azimuth position of the cavity
        PDIS,           // perpendicular distance from symmetric line of cavity gap to machine center
        GAPWIDTH,       // constant gap width of cavity
        PHI0,           // initial phase of cavity
        DESIGNENERGY,   // The mean kinetic energy at exit
        PHASE_MODEL,    // time dependent phase
        AMPLITUDE_MODEL,// time dependent amplitude
        FREQUENCY_MODEL,// time dependent frequency
        SIZE
    };

    /// Exemplar constructor.
    OpalCavity();

    virtual ~OpalCavity();

    /// Make clone.
    virtual OpalCavity *clone(const std::string &name);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalCavity(const OpalCavity &);
    void operator=(const OpalCavity &);

    // Clone constructor.
    OpalCavity(const std::string &name, OpalCavity *parent);

    OpalWake *owk_m;

    BoundaryGeometry *obgeo_m;


};

#endif // OPAL_OpalCavity_HH