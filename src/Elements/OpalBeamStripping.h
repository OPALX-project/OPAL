//
// Class OpalBeamStripping
//   The class of OPAL beam stripping.
//
// Copyright (c) 2018-2019, Pedro Calvo, CIEMAT, Spain
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
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
#ifndef OPAL_OpalBeamStripping_HH
#define OPAL_OpalBeamStripping_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;
// Class OpalBeamStripping
// ------------------------------------------------------------------------
/// The BEAMSTRIPPING element.

class OpalBeamStripping: public OpalElement {

public:

    /// The attributes of class OpalBeamStripping.
    enum {
        PRESSURE = COMMON,   // Pressure in mbar
        TEMPERATURE,         // Temperature of residual gas
        PMAPFN,              // The filename of the mid-plane pressure map
        PSCALE,              // A scalar to scale the P-field
        GAS,                 // Type of gas for residual vaccum
        STOP,                // whether the secondary particles are tracked
        SIZE
    };

    /// Exemplar constructor.
    OpalBeamStripping();

    virtual ~OpalBeamStripping();

    /// Make clone.
    virtual OpalBeamStripping *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC beam stripping.
    virtual void update();

private:

    // Not implemented.
    OpalBeamStripping(const OpalBeamStripping &);
    void operator=(const OpalBeamStripping &);

    // Clone constructor.
    OpalBeamStripping(const std::string &name, OpalBeamStripping *parent);
    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalBeamStripping_HH