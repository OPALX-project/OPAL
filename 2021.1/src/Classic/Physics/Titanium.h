//
// Class Titanium
//   Defines Titanium as a material for particle-matter interactions
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
#ifndef TITANIUM_H
#define TITANIUM_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Titanium: public Material {
    public:
        Titanium():
            Material(22,
                     47.867,
                     4.540,
                     16.16,
                     233.0,
                     {{4.858, 5.489, 5.260e3, 6.511e2, 8.930e-3,
                       4.71, 0.5087, 65.28, 8.806, 5.948}})
        { }
    };
}
#endif
