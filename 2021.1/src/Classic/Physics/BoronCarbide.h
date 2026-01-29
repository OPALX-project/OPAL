//
// Class BoronCarbide
//   Defines BoronCarbide as a material for particle-matter interactions
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
#ifndef BORONCARBIDE_H
#define BORONCARBIDE_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class BoronCarbide: public Material {
    public:
        BoronCarbide():
            Material(26,
                     55.251,
                     2.52,
                     50.13,
                     84.7,
                     {{3.519, 3.963, 6065.0, 1243.0, 7.782e-3,
                       5.013, 0.4707, 85.8, 16.55, 3.211}})
        { }
    };
}
#endif
