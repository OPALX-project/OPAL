/*
 *  Copyright (c) 2025, Jon Thompson
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
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
