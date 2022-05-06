//
// Class GraphiteR6710
//   Defines GraphiteR6710 as a material for particle-matter interactions
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
#ifndef GRAPHITER6710_H
#define GRAPHITER6710_H

#include "Physics/Material.h"

namespace Physics {
    class GraphiteR6710: public Material {
    public:
        GraphiteR6710():
            Material(6,
                     12.0107,
                     1.88,
                     42.70,
                     78.0,
                     {{0.0, 2.601, 1.701e3, 1.279e3, 1.638e-2,
                       3.80133, 0.41590, 12.9966, 117.83, 242.28}})
        { }
    };
}
#endif
