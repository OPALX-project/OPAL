//
// Class Mylar
//   Defines Mylar as a material for particle-matter interactions
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
#ifndef MYLAR_H
#define MYLAR_H

#include "Physics/Material.h"

namespace Physics {
    class Mylar: public Material {
    public:
        Mylar():
            Material(6.702,
                     12.88,
                     1.400,
                     39.95,
                     78.7,
                     {{2.954, 3.350, 1683, 1900, 2.513e-02,
                       1.9259, 0.5550, 27.15125, 26.0665, 6.2768}})
        { }
    };
}
#endif
