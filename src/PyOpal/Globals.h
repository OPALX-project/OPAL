#include <boost/python.hpp>

#include "opal.h"

#ifndef PYOPAL_GLOBALS_H
#define PYOPAL_GLOBALS_H

namespace {
// I get a gcc warning "declared ‘static’ but never defined" if I uncomment this
// but it is not declared static! confused %)
/* 
   void errorHandlerGSL(const char *reason,
                         const char *file,
                         int line,
                         int gsl_errno);
*/
}

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
#endif

#ifdef PYOPAL_GLOBALS_C
    Ippl *ippl = NULL;
    Inform *gmsg = NULL;
#endif
namespace PyOpal {

namespace py = boost::python;

namespace Globals {

void Initialise();
}
}

#endif