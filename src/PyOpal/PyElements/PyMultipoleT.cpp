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
#include "AbsBeamline/MultipoleT.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "Elements/OpalMultipoleT.h"

namespace PyOpal {

    /** Overload this method on the PyOpalObject - it is supposed to be called every
     *  time the cavity is updated (i.e. every time get_field_value is called from
     *  python).
     *
     *  The point is the TimeDependence could have changed and we redo setup just in
     *  case (terribly inefficient, but only effects python side)
     */
    template <>
    void PyOpalObjectNS::PyOpalObject<OpalMultipoleT>::doSetup() {
        object_m->update();
        auto* p = dynamic_cast<MultipoleT*>(object_m->getElement());
        p->initialiseTimeDepencencies();
    }

    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalMultipoleT>::attributes = {
            {"TP", "t_p", "", FLOAT_LIST},
            {"LFRINGE", "left_fringe", "", DOUBLE},
            {"RFRINGE", "right_fringe", "", DOUBLE},
            {"HAPERT", "horizontal_aperture", "", DOUBLE},
            {"VAPERT", "vertical_aperture", "", DOUBLE},
            {"ANGLE", "angle", "", DOUBLE},
            {"EANGLE", "entrance_angle", "", DOUBLE},
            {"MAXFORDER", "maximum_f_order", "", DOUBLE},
            {"MAXXORDER", "maximum_x_order", "", DOUBLE},
            {"ROTATION", "rotation", "", DOUBLE},
            {"VARRADIUS", "variable_radius", "", BOOL},
            {"BBLENGTH", "bounding_box_length", "", DOUBLE},
            {"ENTRYOFFSET", "entry_offset", "", DOUBLE},
            {"L", "length", "", DOUBLE},
            {"DELETEONTRANSVERSEEXIT", "delete_on_transverse_exit", "", BOOL},
            {"SCALING_MODEL", "scaling_model", "", STRING},
        };

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalMultipoleT>::classDocstring =
        "MultipoleT class is a field element that models a Multipole with maxwellian\n"
        "fringe fields. Placement is from the magnet entry point in line with other"
        "elements, however, when VARRADIUS is true this can be switch to about the centre "
        "of the field object by setting ENTRYOFFSET to L/2.\n";

    namespace PyMultipoleT {

        const auto* module_docstring = "multipolet contains the MultipoleT class\n";

        BOOST_PYTHON_MODULE(multipolet) {
            Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalMultipoleT> element;
            auto elementClass = element.make_element_class("MultipoleT");
            element.addGetFieldValue(elementClass, 1.0, 1.0, 1.0, 1e-1);
            elementClass.def("update_time_dependence",
                &PyOpalObjectNS::PyOpalObject<OpalMultipoleT>::doSetup);
        }
    }  // namespace PyMultipoleT
}  // namespace PyOpal
