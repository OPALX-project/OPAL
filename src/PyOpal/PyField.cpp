#include <Python.h>
#include <structmember.h>

#include <boost/python.hpp>

#include "Utilities/OpalException.h"
#include "AbsBeamline/Ring.h"
#include "PyOpal/ExceptionTranslation.h"
#include "PyOpal/PyField.h"

namespace PyOpal {
namespace Field {

py::object get_field_value(double x, double y, double z, double t) {
    std::cerr << "Get Field Value 1" << std::endl;
    Ring* ring = const_cast<Ring*>(Ring::getLastLockedRing());
    std::cerr << "Get Field Value 2" << std::endl;
    if (ring == NULL) {
        std::string err = "Could not find a ring object - maybe a "
           "RingDefinition was not defined or KeepAlive was False";
        throw(OpalException("PyField::get_field_value", err));
    }
    Vector_t R(x, y, z);
    Vector_t P(0, 0, 0);
    Vector_t E, B;
    std::cerr << "Get Field Value 3 " << ring << std::endl;
    int outOfBounds = ring->apply(R, P, t, E, B);
    std::cerr << "Get Field Value 4" << std::endl;
    boost::python::tuple value = boost::python::make_tuple(outOfBounds,
                                          B[0]/10., B[1]/10., B[2]/10.,
                                          E[0], E[1], E[2]);
    std::cerr << "Get Field Value 5" << std::endl;
    return value;
}

BOOST_PYTHON_MODULE(field) {
    ExceptionTranslation::registerExceptions();
    py::scope().attr("__doc__") = field_docstring.c_str();
    py::def("get_field_value",
            get_field_value,
            py::args("x", "y", "z", "t"),
            get_field_value_docstring.c_str()
    );
}

}
}

