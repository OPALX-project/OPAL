//
// Class OpalVacuum
//   Defines the VACUUM element and its attributes.
//
// Copyright (c) 2018 - 2021, Pedro Calvo, CIEMAT, Spain
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
#ifndef OPAL_OpalVacuum_HH
#define OPAL_OpalVacuum_HH

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;

class OpalVacuum: public OpalElement {

public:

    /// The attributes of class OpalVacuum.
    enum {
        GAS = COMMON, // Type of gas for residual vacuum
        PRESSURE,     // Pressure in mbar
        TEMPERATURE,  // Temperature of residual gas
        PMAPFN,       // The filename of the mid-plane pressure map
        PSCALE,       // A scalar to scale the P-field
        STOP,         // whether the secondary particles are tracked
        SIZE
    };

    /// Exemplar constructor.
    OpalVacuum();

    virtual ~OpalVacuum();

    /// Make clone.
    virtual OpalVacuum* clone(const std::string& name);

    /// Update the embedded CLASSIC vacuum.
    virtual void update();

private:

    // Not implemented.
    OpalVacuum(const OpalVacuum&);
    void operator=(const OpalVacuum&);

    // Clone constructor.
    OpalVacuum(const std::string& name, OpalVacuum* parent);

    ParticleMatterInteraction* parmatint_m;
};

#endif // OPAL_OpalVacuum_HH
