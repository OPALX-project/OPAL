//
// Python API for OpalMultipoleT
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

#include "Elements/OpalMultipoleT.h"

namespace PyOpal {

    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalMultipoleT>::attributes = {
            {"TP", "t_p", "", PyOpalObjectNS::FLOAT_LIST},
            {"LFRINGE", "left_fringe", "", PyOpalObjectNS::DOUBLE},
            {"RFRINGE", "right_fringe", "", PyOpalObjectNS::DOUBLE},
            {"HAPERT", "horizontal_aperture", "", PyOpalObjectNS::DOUBLE},
            {"VAPERT", "vertical_aperture", "", PyOpalObjectNS::DOUBLE},
            {"ANGLE", "angle", "", PyOpalObjectNS::DOUBLE},
            {"EANGLE", "entrance_angle", "", PyOpalObjectNS::DOUBLE},
            {"MAXFORDER", "maximum_f_order", "", PyOpalObjectNS::DOUBLE},
            {"MAXXORDER", "maximum_x_order", "", PyOpalObjectNS::DOUBLE},
            {"ROTATION", "rotation", "", PyOpalObjectNS::DOUBLE},
            {"VARRADIUS", "variable_radius", "", PyOpalObjectNS::BOOL},
            {"BBLENGTH", "bounding_box_length", "", PyOpalObjectNS::DOUBLE},
            {"ENTRYOFFSET", "entry_offset", "", PyOpalObjectNS::DOUBLE},
            {"L", "length", "", PyOpalObjectNS::DOUBLE}};

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalMultipoleT>::classDocstring =
        "MultipoleT class is a field element that models a Multipole with maxwellian\n"
        "fringe fields. Placement is from the magnet entry point in line with other"
        "elements, however, when VARRADIUS is true this can be switch to about the centre "
        "of the field object by setting ENTRYOFFSET to L/2.\n";

    namespace PyMultipoleT {

        const char* module_docstring = "multipolet contains the MultipoleT class\n";

        BOOST_PYTHON_MODULE(multipolet) {
            PyOpal::Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalMultipoleT> element;
            auto elementClass = element.make_element_class("MultipoleT");
            element.addGetFieldValue(elementClass, 1.0, 1.0, 1.0, 1e-1);
        }
    }  // namespace PyMultipoleT
}  // namespace PyOpal
