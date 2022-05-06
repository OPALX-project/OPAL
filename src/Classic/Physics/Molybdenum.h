//
// Class Molybdenum
//   Defines Molybdenum as a material for particle-matter interactions
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
#ifndef MOLYBDENUM_H
#define MOLYBDENUM_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Molybdenum: public Material {
    public:
        Molybdenum():
            Material(42,
                     95.95,
                     10.22,
                     9.80,
                     424.0,
                     {{6.424, 7.248, 9.545e3, 4.802e2, 5.376e-3,
                       9.276, 0.418, 157.1, 8.038, 1.29}})
        { }
    };
}
#endif
