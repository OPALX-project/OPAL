#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "AbsBeamline/Ring.h" // classic
#include "Utilities/RingSection.h" // classic
#include "Physics/Units.h" // classic
#include "Elements/OpalRingDefinition.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyRingDefinition {

/** RingDefinition handles access to the Ring object.
 * 
 *  Placements in the Ring object are handled by RingSection. To keep the UI
 *  simple I don't give user direct access to RingSection, but provide a few
 *  convenience methods here. 
 * 
 *  General note that I am following the behaviour of Opal-T and
 *  Opal-Cyclotron which go through a "line" object for element placements. But
 *  the "Line" doesn't know about the actual position of the elements, so poor
 *  user has to do some kludge where they do field lookups and lookup position
 *  of elements from Ring, which can only be done after the "Line" setup is 
 *  done (which in itself is done by Track). It's a bit intricate and arcane, 
 *  would be better to access directly through Ring, but that's not what we've
 *  done.
 */

std::string ring_definition_docstring = std::string();

const char* module_docstring = 
"ring_definition module contains the RingDefinition class";

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalRingDefinition>::classDocstring =
"RingDefinition class enables definition of the global ring parameters\n"
"\n"
"RingDefinition holds definitions of the global ring parameters like the\n"
"position of the first lattice element and position of the beam.\n"
"Internally, the RingDefinition handles placement of OPAL objects. Once they\n"
"have been placed, (e.g. Track has been executed), it is possible to\n"
"interrogate the RingDefinition to find what are the element positions as\n"
"placed by PyOpal, for debugging purposes.\n"
"\n"
"See also pyopal.objects.line.Line class, which handles the sequential\n"
"element placements.\n\n";

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

BOOST_PYTHON_MODULE(ring_definition) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalRingDefinition> element;
    auto elementClass = element.make_element_class("RingDefinition");
    element.addGetFieldValue(elementClass, 1.0, 1.0, 1.0, 1e-1);
}

}
}
