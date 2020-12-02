#include "Classic/Utilities/ClassicException.h"
#include "PyOpal/ExceptionTranslation.h"

namespace PyOpal {
namespace ExceptionTranslation {
void registerExceptions() {
    // handle std::exception (e.g. for IO errors)
    py::register_exception_translator<std::exception>(translateException<std::exception>);
    // handle Opal exceptions (they all inherit from ClassicException)
    py::register_exception_translator<ClassicException>(translateOpalException<ClassicException>);
}
}
}