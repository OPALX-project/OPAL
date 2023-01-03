#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Elements/OpalRingDefinition.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyRingDefinition {

std::string ring_definition_docstring = std::string();


const char* module_docstring = "build a ring_definition";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalRingDefinition>::attributes = {
    {"LAT_RINIT", "lattice_initial_r", "", PyOpalObjectNS::DOUBLE},
    {"LAT_PHIINIT", "lattice_initial_phi", "", PyOpalObjectNS::DOUBLE},
    {"LAT_THETAINIT", "lattice_initial_theta", "", PyOpalObjectNS::DOUBLE},
    {"BEAM_RINIT", "beam_initial_r", "", PyOpalObjectNS::DOUBLE},
    {"BEAM_PHIINIT", "beam_initial_phi", "", PyOpalObjectNS::DOUBLE},
    {"BEAM_PRINIT", "beam_initial_pr", "", PyOpalObjectNS::DOUBLE},
    {"HARMONIC_NUMBER", "harmonic_number", "", PyOpalObjectNS::DOUBLE},
    {"SYMMETRY", "symmetry", "", PyOpalObjectNS::INT},
    {"SCALE", "scale", "", PyOpalObjectNS::DOUBLE},
    {"RFFREQ", "rf_frequency", "", PyOpalObjectNS::DOUBLE},
    {"IS_CLOSED", "is_closed", "", PyOpalObjectNS::BOOL}, // BUG in underlying code
    {"MIN_R", "minimum_r", "", PyOpalObjectNS::DOUBLE},
    {"MAX_R", "maximum_r", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalRingDefinition>::classDocstring = "";

BOOST_PYTHON_MODULE(ring_definition) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalRingDefinition> element;
    auto elementClass = element.make_element_class("RingDefinition");
    element.addGetFieldValue(elementClass, 1.0, 1.0, 1.0, 1e-1);
}

}
}
