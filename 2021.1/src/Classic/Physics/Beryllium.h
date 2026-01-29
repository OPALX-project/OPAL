//
// Class Beryllium
//   Defines Beryllium as a material for particle-matter interactions
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
#ifndef BERYLLIUM_H
#define BERYLLIUM_H

#include "Physics/Material.h"

namespace Physics {
    class Beryllium: public Material {
    public:
        Beryllium():
            Material(4,
                     9.0121831,
                     1.848,
                     65.19,
                     63.7,
                     {{2.248, 2.590, 9.660e2, 1.538e2, 3.475e-2,
                       2.1895, 0.47183, 7.2362, 134.30, 197.96}})
        { }
    };
}
#endif
