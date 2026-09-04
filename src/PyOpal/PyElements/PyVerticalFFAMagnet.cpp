//
// Python API for OpalVerticalFFAMagnet
//
// Copyright (c) 2022 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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

#include "Elements/OpalVerticalFFAMagnet.h"

namespace PyOpal {
    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet>::attributes = {
            {"B0", "b0", "", PyOpalObjectNS::DOUBLE},
            {"FIELD_INDEX", "field_index", "", PyOpalObjectNS::DOUBLE},
            {"WIDTH", "width", "", PyOpalObjectNS::DOUBLE},
            {"MAX_HORIZONTAL_POWER", "max_horizontal_power", "", PyOpalObjectNS::INT},
            {"END_LENGTH", "end_length", "", PyOpalObjectNS::DOUBLE},
            {"CENTRE_LENGTH", "centre_length", "", PyOpalObjectNS::DOUBLE},
            {"BB_LENGTH", "bb_length", "", PyOpalObjectNS::DOUBLE},
            {"HEIGHT_NEG_EXTENT", "height_neg_extent", "", PyOpalObjectNS::DOUBLE},
            {"HEIGHT_POS_EXTENT", "height_pos_extent", "", PyOpalObjectNS::DOUBLE},
    };

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet>::classDocstring =
        "VerticalFFAMagnet class is a field element that models a Vertical FFA magnet.";

    namespace PyVerticalFFAMagnet {

        const char* module_docstring = "vertical_ffa_magnet contains the VerticalFFAMagnet class";

        BOOST_PYTHON_MODULE(vertical_ffa_magnet) {
            PyOpal::Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet> element;
            auto elementClass = element.make_element_class("VerticalFFAMagnet");
            element.addGetFieldValue(elementClass, 1.0, 1.0, 1.0, 1e-1);
        }

    }  // namespace PyVerticalFFAMagnet
}  // namespace PyOpal
