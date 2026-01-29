//
// Class AluminaAL2O3
//   Defines Alumina as a material for particle-matter interactions
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
#ifndef ALUMINAAL2O3_H
#define ALUMINAAL2O3_H

#include "Physics/Material.h"

#include <cmath>

namespace Physics {
    class AluminaAL2O3: public Material {
    public:
        AluminaAL2O3():
            Material(50,
                     101.9600768,
                     3.97,
                     27.94,
                     145.2,
                     {{1.187e1, 1.343e1, 1.069e4, 7.723e2, 2.153e-2,
                       5.4, 0.53, 103.1, 3.931, 7.767}})
        { }
    };
}
#endif
