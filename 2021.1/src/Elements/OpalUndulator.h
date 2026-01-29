//
// Class OpalUndulator
//   Defines the Undulator/Wiggler element and its attributes.
//
// Copyright (c) 2020, Arnau Albà, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// Implemented as part of the MSc thesis
// "Start-to-End Modelling of the AWA Micro-Bunched Electron Cooling POP-Experiment"
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
#ifndef OPAL_OpalUndulator_HH
#define OPAL_OpalUndulator_HH

#include "Elements/OpalElement.h"

class OpalUndulator : public OpalElement {
public:
    enum {
        K = COMMON,      // The undulator parameter
        LAMBDA,          // The undulator period
        NUMPERIODS,      // Number of undulator periods
        ANGLE,           // Polarisation angle of the undulator
        FNAME,           // File specifying the wanted output from the full wave simulation
        MESHLENGTH,      // Size of the computational domain
        MESHRESOLUTION,  // Size of the grid-cells
        TRUNORDER,       // Order of the Absorbing Boundary Conditions, 1st or 2nd
        TOTALTIME,       // Total time of the full wave simulation
        DTBUNCH,         // Time-step for particle update can be smaller than field update step
        SIZE
    };

    /// Exemplar constructor.
    OpalUndulator();

    virtual ~OpalUndulator();

    /// Make clone.
    virtual OpalUndulator* clone(const std::string& name);

    /// Update the embedded CLASSIC drift.
    virtual void update();

private:
    // Not implemented.
    OpalUndulator(const OpalUndulator&);
    void operator=(const OpalUndulator&);

    // Clone constructor.
    OpalUndulator(const std::string& name, OpalUndulator* parent);
};

#endif  // OPAL_OpalUndulator_HH
