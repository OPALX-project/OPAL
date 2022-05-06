//
// Class Copper
//   Defines Copper as a material for particle-matter interactions
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
#ifndef COPPER_H
#define COPPER_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Copper: public Material {
    public:
        Copper():
            Material(29,
                     63.546,
                     8.96,
                     12.86,
                     322.0,
                     {{3.969, 4.194, 4.649e3, 8.113e1, 2.242e-2,
                       3.114, 0.5236, 76.67, 7.62, 6.385}})
        { }
    };
}
#endif
