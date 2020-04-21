//
// Header file AmrMultiGridCore
//   Includes all AMR solver core headers.
//
// Copyright (c) 2017 - 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#ifndef AMR_MULTI_GRID_CORE_H
#define AMR_MULTI_GRID_CORE_H

// boundary handlers
#include "AmrDirichletBoundary.h"
#include "AmrOpenBoundary.h"
#include "AmrPeriodicBoundary.h"

// interpolaters
#include "AmrTrilinearInterpolater.h"
#include "AmrLagrangeInterpolater.h"
#include "AmrPCInterpolater.h"

// base level solvers
#include "BottomSolver.h"
#include "BelosBottomSolver.h"
#include "Amesos2BottomSolver.h"
#include "MueLuBottomSolver.h"

#include "AmrSmoother.h"

// preconditioners
#include "Ifpack2Preconditioner.h"
#include "MueLuPreconditioner.h"

#endif
