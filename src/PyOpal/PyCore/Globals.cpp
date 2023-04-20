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

#include <cstring>
#include <gsl/gsl_errno.h>
#include "Utilities/OpalException.h"

// Note the gymnastics here - we only want to define gmsg and ippl once
#define PYOPAL_GLOBALS_C
#include "PyOpal/PyCore/Globals.h"

namespace {
    void errorHandlerGSL(const char *reason,
                         const char *file,
                         int line,
                         int gsl_errno) {
        throw OpalException(file, reason);
        if (line || gsl_errno) {;} // disable gcc warning; does nothing
    }
}

namespace PyOpal {
namespace Globals {
void Initialise() {
    if (gmsg == nullptr) {
        gmsg = new Inform("OPAL");
    }
    if (ippl == nullptr) {
        int argc = 3;
        // I am not strong on the C-style strings, but if I understand correctly
        // there is a secret null pointer at the end of each string, hence the
        // char arrays have to be one character longer than you might think.
        char* argvr[4];
        argvr[0] = new char[7];
        strcpy(argvr[0], "pyopal");
        argvr[1] = new char[12];
        strcpy(argvr[1], "--processes");
        argvr[2] = new char[2];
        strcpy(argvr[2], "3");
        argvr[3] = nullptr;
        char** argv = argvr;
        // Ippl is a typedef of IpplInfo in ippl/Utilities
        ippl = new Ippl(argc, argv);
    }
    gsl_set_error_handler(&errorHandlerGSL);
}
}
}
