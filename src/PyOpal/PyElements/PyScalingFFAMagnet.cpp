//
// Python API for OpalScalingFFAMagnet
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

#include "AbsBeamline/ScalingFFAMagnet.h"
#include "Elements/OpalScalingFFAMagnet.h"

namespace PyOpal {

    /** Overload this method on the PyOpalObject - it is supposed to be called every
     *  time the magnet is updated (i.e. every time get_field_value is called from python).
     *
     *  The point is the EndFieldModel could have changed and we redo setup just in case
     *  (terribly inefficient)
     */
    template <>
    void PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet>::doSetup() {
        ScalingFFAMagnet* magnet = dynamic_cast<ScalingFFAMagnet*>(object_m->getElement());
        magnet->setupEndField();
    }

    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet>::attributes = {
            {"B0", "b0", "", PyOpalObjectNS::DOUBLE},
            {"R0", "r0", "", PyOpalObjectNS::DOUBLE},
            {"FIELD_INDEX", "field_index", "", PyOpalObjectNS::DOUBLE},
            {"TAN_DELTA", "tan_delta", "", PyOpalObjectNS::DOUBLE},
            {"MAX_Y_POWER", "max_vertical_power", "", PyOpalObjectNS::INT},
            {"END_FIELD_MODEL", "end_field_model", "", PyOpalObjectNS::STRING},
            {"END_LENGTH", "end_length", "", PyOpalObjectNS::DOUBLE},
            {"CENTRE_LENGTH", "centre_length", "", PyOpalObjectNS::DOUBLE},
            {"HEIGHT", "height", "", PyOpalObjectNS::DOUBLE},
            {"RADIAL_NEG_EXTENT", "radial_neg_extent", "", PyOpalObjectNS::DOUBLE},
            {"RADIAL_POS_EXTENT", "radial_pos_extent", "", PyOpalObjectNS::DOUBLE},
            {"MAGNET_START", "magnet_start", "", PyOpalObjectNS::DOUBLE},
            {"MAGNET_END", "magnet_end", "", PyOpalObjectNS::DOUBLE},
            {"AZIMUTHAL_EXTENT", "azimuthal_extent", "", PyOpalObjectNS::DOUBLE},
    };

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet>::classDocstring =
        "ScalingFFAMagnet class is a field element that models a Scaling FFA magnet.";

    namespace PyScalingFFAMagnet {

        const char* module_docstring = "scaling_ffa_magnet contains the ScalingFFAMagnet class";

        std::string update_docstr =
            "Check for changes to the EndFieldModel and update the ScalingFFAMagnet "
            "appropriately.\n"
            "This is done automatically the first time the ScalingFFAMagnet is used but not for\n"
            "subsequent uses. WARNING: if user changes the end field model, user must call\n"
            "'update_end_field' manually to load the new parameters.\n"
            "\n"
            "May throw RuntimeError if the EndFieldModel is not valid or cannot be found.\n"
            "\n"
            "Returns None.\n";

        BOOST_PYTHON_MODULE(scaling_ffa_magnet) {
            PyOpal::Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet> element;
            auto elementClass = element.make_element_class("ScalingFFAMagnet");
            element.addGetFieldValue(elementClass, 1.0, 1.0, 1.0, 1e-1);

            elementClass.def(
                "update_end_field", &PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet>::doSetup,
                update_docstr.c_str());
        }

    }  // namespace PyScalingFFAMagnet
}  // namespace PyOpal
