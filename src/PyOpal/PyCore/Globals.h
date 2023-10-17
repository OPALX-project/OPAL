//
// Python API for PolynomialCoefficient (part of the multidimensional polynomial fitting routines)
//
// Copyright (c) 2008-2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//

#ifndef PYOPAL_PYCORE_GLOBALS_H
#define PYOPAL_PYCORE_GLOBALS_H
#include "opal.h"

// Ippl and gmsg should only be built once, in globals.cc.o
//
// We declare them as extern for all other code to be consistent with usage in 
// C++ Opal; but we define them if we are including from Globals.cc so that the 
// extern is satisfied.
//
// Maybe this implements one gmsg per python module - which is wrong? So maybe
// some more linker/cmake dark arts required here? Do they have much global 
// state?
#ifndef PYOPAL_GLOBALS_C
    extern Ippl *ippl;
    extern Inform *gmsg;
    extern Inform *gmsgALL;
#endif

#ifdef PYOPAL_GLOBALS_C
    Ippl *ippl = nullptr;
    Inform *gmsg = nullptr;
    Inform *gmsgALL = nullptr;
#endif
namespace PyOpal {
namespace Globals {
/** Globals namespace provides routines to initialise global objects:
 *  - ippl
 *  - gmsg
 */
void Initialise();
}
}

#endif
