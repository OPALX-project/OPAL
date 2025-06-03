//
// Python API for OpalVariableRFCavity
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

#include "AbsBeamline/VariableRFCavity.h"
#include "Elements/OpalVariableRFCavity.h"
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
    void PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity>::doSetup() {
        object_m->update();
        VariableRFCavity* rf = dynamic_cast<VariableRFCavity*>(object_m->getElement());
        rf->initialise();
    }

    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity>::attributes = {
            {"PHASE_MODEL", "phase_model", "", PyOpalObjectNS::STRING},
            {"AMPLITUDE_MODEL", "amplitude_model", "", PyOpalObjectNS::STRING},
            {"FREQUENCY_MODEL", "frequency_model", "", PyOpalObjectNS::STRING},
            {"WIDTH", "width", "", PyOpalObjectNS::DOUBLE},
            {"HEIGHT", "height", "", PyOpalObjectNS::DOUBLE},
            {"L", "length", "", PyOpalObjectNS::DOUBLE},
    };

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity>::classDocstring =
        "VariableRFCavity class is a field element that models a rectangular RF cavity.";

    namespace PyVariableRFCavity {

        const char* module_docstring = "variable_rf_cavity contains the VariableRFCavity class";

        BOOST_PYTHON_MODULE(variable_rf_cavity) {
            PyOpal::Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity> element;
            auto elementClass = element.make_element_class("VariableRFCavity");
            element.addGetFieldValue(elementClass, 1.0, Units::s2ns, 1.0, 1e-1);
            elementClass.def(
                "update_time_dependence",
                &PyOpalObjectNS::PyOpalObject<OpalVariableRFCavity>::doSetup);
        }
    }  // namespace PyVariableRFCavity
}  // namespace PyOpal
