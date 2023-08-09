#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Classic/AbsBeamline/VariableRFCavityFringeField.h"
#include "Elements/OpalVariableRFCavityFringeField.h"

namespace PyOpal {

/** Overload this method on the PyOpalObject - it is supposed to be called every
 *  time the cavity is updated (i.e. every time get_field_value is called from 
 *  python).
 * 
 *  The point is the TimeDependence could have changed and we redo setup just in
 *  case (terribly inefficient, but only effects python side)
 */
template <>
void PyOpalObjectNS::PyOpalObject<OpalVariableRFCavityFringeField>::doSetup() {
    object_m->update();
    VariableRFCavityFringeField* rf =
                      dynamic_cast<VariableRFCavityFringeField*>(object_m->getElement());
    rf->initialise();
}

namespace PyVariableRFCavityFringeField {

const char* module_docstring = 
    "variable_rf_cavity contains the VariableRFCavityFringeField class";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalVariableRFCavityFringeField>::attributes = {
    {"PHASE_MODEL", "phase_model", "", PyOpalObjectNS::STRING},
    {"AMPLITUDE_MODEL", "amplitude_model", "", PyOpalObjectNS::STRING},
    {"FREQUENCY_MODEL", "frequency_model", "", PyOpalObjectNS::STRING},
    {"WIDTH", "width", "", PyOpalObjectNS::DOUBLE},
    {"HEIGHT", "height", "", PyOpalObjectNS::DOUBLE},
    {"CENTRE_LENGTH", "centre_length", "", PyOpalObjectNS::DOUBLE},
    {"END_LENGTH", "end_length", "", PyOpalObjectNS::DOUBLE},
    {"CAVITY_CENTRE", "cavity_centre", "", PyOpalObjectNS::DOUBLE},
    {"MAX_ORDER", "max_order", "", PyOpalObjectNS::DOUBLE},
    {"L", "length", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalVariableRFCavityFringeField>::classDocstring =
"VariableRFCavityFringeField class is a field element that models a rectangular RF cavity.";

BOOST_PYTHON_MODULE(variable_rf_cavity_fringe_field) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalVariableRFCavityFringeField> element;
    auto elementClass = element.make_element_class("VariableRFCavityFringeField");
    element.addGetFieldValue(elementClass, 1e+3, 1e+9, 1.0, 1e-1);
    elementClass.def("update_time_dependence",
                &PyOpalObjectNS::PyOpalObject<OpalVariableRFCavityFringeField>::doSetup);
}
}
}
