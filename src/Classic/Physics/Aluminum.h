//
// Class Aluminum
//   Defines Aluminum as a material for particle-matter interactions
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
#ifndef ALUMINUM_H
#define ALUMINUM_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class Aluminum: public Material {
    public:
        Aluminum():
            Material(13,
                     26.9815384,
                     2.699,
                     24.01,
                     166.0,
                     {{4.154, 4.739, 2.766e3, 1.645e2, 2.023e-2,
                       2.5, 0.625, 45.7, 0.1, 4.359}})
        { }
    };
}
#endif
