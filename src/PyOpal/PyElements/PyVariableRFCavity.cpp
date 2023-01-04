#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Classic/AbsBeamline/VariableRFCavity.h"
#include "Elements/OpalVariableRFCavity.h"

namespace PyOpal {

/** Overload this method on the PyOpalObject - it is supposed to be called every
 *  time the cavity is updated (i.e. every time get_field_value is called from 
 *  python).
 * 
 *  The point is the TimeDependence could have changed and we redo setup just in
 *  case (terribly inefficient, but only effects python side)
 */
template <>
void PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity>::doSetup() {
    VariableRFCavity* rf =
                      dynamic_cast<VariableRFCavity*>(object_m->getElement());
    rf->initialise();
}

namespace PyVariableRFCavity {

const char* module_docstring = 
    "variable_rf_cavity contains the VariableRFCavity class";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity>::attributes = {
    {"PHASE_MODEL", "phase_model", "", PyOpalObjectNS::STRING},
    {"AMPLITUDE_MODEL", "amplitude_model", "", PyOpalObjectNS::STRING},
    {"FREQUENCY_MODEL", "frequency_model", "", PyOpalObjectNS::STRING},
    {"WIDTH", "variable_radius", "", PyOpalObjectNS::DOUBLE},
    {"HEIGHT", "bounding_box_length", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity>::classDocstring = 
"VariableRFCavity class is a field element that models a rectangular RF cavity.";

BOOST_PYTHON_MODULE(variable_rf_cavity) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity> element;
    auto elementClass = element.make_element_class("VariableRFCavity");
    element.addGetFieldValue(elementClass, 1e+3, 1.0, 1.0, 1e-1);
}
}
}
