//
// File Steppers
//   Integrator enum.
//
// Copyright (c) 2017, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef STEPPERS_H
#define STEPPERS_H

#include "Steppers/RK4.h"
#include "Steppers/LF2.h"

namespace Steppers {
    enum TimeIntegrator: short {
        UNDEFINED   = -1,
        RK4         = 0,
        LF2         = 1,
        MTS         = 2
    };
};

#endif
