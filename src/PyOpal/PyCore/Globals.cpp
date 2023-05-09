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
#include <Python.h>
#include <boost/python.hpp>

#include "Utility/IpplInfo.h" // that is ippl land
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
void printArgv(char** argv) { // debugging
    int i = 0;
    while (argv[i] != nullptr) {
        std::cerr << i << " " << std::string(argv[i]) << std::endl;
        ++i;
    }

}

void Initialise() {
    if (ippl == nullptr) {
        PyObject* pyargv = PySys_GetObject("argv"); // this is a borrowed ref
        boost::python::handle<> wrapper(boost::python::borrowed(pyargv)); // now wrapper owns the ref
        boost::python::list myList(wrapper);

        int argc = int(boost::python::len(myList));
        // I am not strong on the C-style strings, but if I understand correctly
        // there is a secret null pointer at the end of each string, hence the
        // char arrays have to be one character longer than you might think.
        char* argvr[argc+1];
        argvr[0] = new char[7];
        strcpy(argvr[0], "pyopal");
        for (int i = 1; i < argc; ++i) {
            int stringLength(boost::python::len(myList[i]));
            const char* value = boost::python::extract<const char*>(
                                                 boost::python::str(myList[i]));
            argvr[i] = new char[stringLength+1];
            strcpy(argvr[i], value);
        }
        argvr[argc] = nullptr;
        // and here is another secret nullptr to mark the end of the array of strings
        // don't forget it, you might spend days debugging the segv otherwise...
        char** argv = argvr;
        // Ippl is a typedef of IpplInfo in ippl/Utilities
        ippl = new Ippl(argc, argv);
    }
    if (gmsg == nullptr) {
        IpplInfo::instantiateGlobals();
        gmsg = new Inform("OPAL");
        gmsgALL = new Inform("OPAL", INFORM_ALL_NODES);
    }
    gsl_set_error_handler(&errorHandlerGSL);
}
}
}
