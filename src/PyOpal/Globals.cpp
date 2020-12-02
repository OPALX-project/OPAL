#include <gsl/gsl_errno.h>

#include "Utilities/OpalException.h"
#define PYOPAL_GLOBALS_C
#include "PyOpal/Globals.h"

namespace {
    void errorHandlerGSL(const char *reason,
                         const char *file,
                         int /*line*/,
                         int /*gsl_errno*/) {
        throw OpalException(file, reason);
    }
}

namespace PyOpal {
namespace Globals {
void Initialise() {
    if (ippl == NULL) {
        //ippl = new Ippl(argc, argv);
    }
    if (gmsg == NULL) {
        gmsg = new  Inform("OPAL");
    }
    gsl_set_error_handler(&errorHandlerGSL);
}
}
}