// Copyright (c) 2023, Chris Rogers
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "Elements/OpalOutputPlane.h"

namespace PyOpal {
    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalOutputPlane>::attributes = {
            {"CENTRE", "centre", "", PyOpalObjectNS::FLOAT_LIST},
            {"NORMAL", "normal", "", PyOpalObjectNS::FLOAT_LIST},
            {"XSTART", "x_start", "", PyOpalObjectNS::DOUBLE},
            {"XEND", "x_end", "", PyOpalObjectNS::DOUBLE},
            {"YSTART", "y_start", "", PyOpalObjectNS::DOUBLE},
            {"YEND", "y_end", "", PyOpalObjectNS::DOUBLE},
            {"PLACEMENT_STYLE", "placement_style", "", PyOpalObjectNS::PREDEFINED_STRING},
            {"ALGORITHM", "algorithm", "", PyOpalObjectNS::PREDEFINED_STRING},
            {"TOLERANCE", "tolerance", "", PyOpalObjectNS::DOUBLE},
            {"REFERENCE_ALIGNMENT_PARTICLE", "reference_alignment_particle", "",
             PyOpalObjectNS::INT},
            {"OUTFN", "output_filename", "",
             PyOpalObjectNS::STRING},  // OUTFN comes from OpalElement (yes, all elements can have a
                                       // filename!)
            {"VERBOSE", "verbose_level", "", PyOpalObjectNS::INT},
            {"WIDTH", "width", "", PyOpalObjectNS::DOUBLE},
            {"HEIGHT", "height", "", PyOpalObjectNS::DOUBLE},
            {"RADIUS", "radius", "", PyOpalObjectNS::DOUBLE}};

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalOutputPlane>::classDocstring =
        "OutputPlane is used to generate output data based on particle tracks crossing a plane.";

    namespace PyOutputPlane {

        const char* module_docstring = "output_plane contains the OutputPlane class";

        BOOST_PYTHON_MODULE(output_plane) {
            PyOpal::Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalOutputPlane> element;
            auto elementClass = element.make_element_class("OutputPlane");
        }

    }  // namespace PyOutputPlane
}  // namespace PyOpal
