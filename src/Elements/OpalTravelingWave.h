//
// Class OpalTravelingWave
//   The TRAVELINGWAVE element.
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
#ifndef OPAL_OpalTravelingWave_HH
#define OPAL_OpalTravelingWave_HH

#include "Elements/OpalElement.h"

class OpalWake;

class OpalTravelingWave: public OpalElement {

public:

    /// The attributes of class OpalTravelingWave.
    enum {
        VOLT = COMMON,  // The peak voltage.
        DVOLT,          // The peak voltage error
        FREQ,           // The RF frequency.
        LAG,            // The phase lag.
        DLAG,           // The phase lag error
        FMAPFN,         // The filename of the fieldmap
        APVETO,         // Do not use this cavity in the Autophase procedure
        FAST,           // Faster but less accurate
        NUMCELLS,       // Number of cells in a TW structure
        DESIGNENERGY,   // The mean kinetic energy at exit
        MODE,           // The phase shift between cells
        SIZE
    };

    /// Exemplar constructor.
    OpalTravelingWave();

    virtual ~OpalTravelingWave();

    /// Make clone.
    virtual OpalTravelingWave *clone(const std::string &name);

    /// Update the embedded CLASSIC cavity.
    virtual void update();

private:

    // Not implemented.
    OpalTravelingWave(const OpalTravelingWave &);
    void operator=(const OpalTravelingWave &);

    // Clone constructor.
    OpalTravelingWave(const std::string &name, OpalTravelingWave *parent);

    OpalWake *owk_m;
};

#endif // OPAL_OpalTravelingWave_HH