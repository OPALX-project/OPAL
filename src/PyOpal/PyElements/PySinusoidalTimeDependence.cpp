//
// Class PySinusoidalTimeDependence
//   Python interface for the time dependence class that generates sine waves
//
// Copyright (c) 2025, Jon Thompson, Isis, RAL, UK
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

#include "Elements/OpalSinusoidalTimeDependence.h"
#include "PyOpal/PyElements/PyAbstractTimeDependence.h"

namespace PyOpal {
    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalSinusoidalTimeDependence>::attributes = {
            {"FREQUENCIES", "frequencies", "", FLOAT_LIST},
            {"AMPLITUDES", "amplitudes", "", FLOAT_LIST},
            {"PHASE_OFFSETS", "phase_offsets", "", FLOAT_LIST},
            {"DC_OFFSETS", "dc_offsets", "", FLOAT_LIST},
    };

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalSinusoidalTimeDependence>::classDocstring =
        "SinusoidalTimeDependence class enables modelling time dependent behaviour.\n"
        "\n"
        "SinsoidalTimeDependence class enables modelling lattice parameters that are\n"
        "changing in time, for example RF voltages, frequencies, etc. The dependence\n"
        "is referenced and set to field elements by use of the OPAL name (note not\n"
        "any python name). The OPAL name is set/retrieved using 'set_opal_name' and\n"
        "'get_opal_name' respectively.\n";

    namespace PyOpalSinusoidalTimeDependence {

        const char* module_docstring =
            "sinusoidal_time_dependence module holds a SinusoidalTimeDependence class";

        BOOST_PYTHON_MODULE(sinusoidal_time_dependence) {
            Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalSinusoidalTimeDependence> element;
            auto elementClass = element.make_class("SinusoidalTimeDependence");
            elementClass.def(
                "function",
                &PyOpal::PyAbstractTimeDependence::function<OpalSinusoidalTimeDependence>);
            elementClass.def("update", &PyOpalObjectNS::update<OpalSinusoidalTimeDependence>);
        }

    }  // namespace PyOpalSinusoidalTimeDependence
}  // namespace PyOpal
