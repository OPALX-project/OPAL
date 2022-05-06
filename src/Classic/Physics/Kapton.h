//
// Class Kapton
//   Defines Kapton as a material for particle-matter interactions
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
#ifndef KAPTON_H
#define KAPTON_H

#include "Physics/Material.h"

namespace Physics {
    class Kapton: public Material {
    public:
        Kapton():
            Material(6,
                     12,
                     1.420,
                     40.58,
                     79.6,
                     {{0.0, 2.601, 1.701e3, 1.279e3, 1.638e-2,
                       3.83523, 0.42993, 12.6125, 227.41, 188.97}})
        { }
    };
}
#endif
