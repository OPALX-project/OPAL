//
// Class Water
//   Defines Water as a material for particle-matter interactions
//
//   Data taken from Standard Atomic Weight 2019 (https://www.qmul.ac.uk/sbcs/iupac/AtWt/), 
//   atomic properties from PDG database (https://pdg.lbl.gov/2020/AtomicNuclearProperties)
//   and fit coefficients from ICRU Report 49.
//
// Copyright (c) 2019 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef WATER_H
#define WATER_H

#include "Physics/Material.h"

namespace Physics {
    class Water: public Material {
    public:
        Water():
            Material(10,
                     18.0152,
                     1,
                     36.08,
                     75.0,
                     {{4.015, 4.542, 3.955e3, 4.847e2, 7.904e-3,
                       2.9590, 0.53255, 34.247, 60.655, 15.153}})
        { }
    };
}
#endif
