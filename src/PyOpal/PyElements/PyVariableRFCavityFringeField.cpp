//
// Python API for OpalVariableRFCavityFringeField
//
// Copyright (c) 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "Classic/AbsBeamline/VariableRFCavityFringeField.h"
#include "Elements/OpalVariableRFCavityFringeField.h"
#include "Physics/Units.h"

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

    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalVariableRFCavityFringeField>::attributes = {
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

    namespace PyVariableRFCavityFringeField {

        const char* module_docstring =
            "variable_rf_cavity contains the VariableRFCavityFringeField class";

        BOOST_PYTHON_MODULE(variable_rf_cavity_fringe_field) {
            PyOpal::Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalVariableRFCavityFringeField> element;
            auto elementClass = element.make_element_class("VariableRFCavityFringeField");
            element.addGetFieldValue(elementClass, 1.0, Units::s2ns, 1.0, 1e-1);
            elementClass.def(
                "update_time_dependence",
                &PyOpalObjectNS::PyOpalObject<OpalVariableRFCavityFringeField>::doSetup);
        }
    }  // namespace PyVariableRFCavityFringeField
}  // namespace PyOpal
