#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Elements/OpalElement.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyOpalElement {

std::string track_run_docstring = std::string();


const char* module_docstring = "opal element base class";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalElement>::attributes = {
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalElement>::classDocstring = 
"OpalElement is a base class for Opal element objects. In particular it is\n"
"used by Line class to handle wrapping of objects";

template <>
PyOpalObjectNS::PyOpalObject<OpalElement>::PyOpalObject() : object_m(NULL) {}

BOOST_PYTHON_MODULE(opal_element) {
    ExceptionTranslation::registerExceptions();
    PyOpal::Globals::Initialise();
    PyOpalObjectNS::PyOpalObject<OpalElement> element;
    element.make_element_class("OpalElement");
}

}
}