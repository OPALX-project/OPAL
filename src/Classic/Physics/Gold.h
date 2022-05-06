//
// Class Gold
//   Defines Gold as a material for particle-matter interactions
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
#ifndef GOLD_H
#define GOLD_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Gold: public Material {
    public:
        Gold():
            Material(79,
                     196.966570,
                     19.32,
                     6.46,
                     790.0,
                     {{4.844, 5.458, 7.852e3, 9.758e2, 2.077e-2,
                       3.223, 0.5883, 232.7, 2.954, 1.05}})
        { }
    };
}
#endif
