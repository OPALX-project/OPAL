//
// Python API for TrackCmd
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

#include "Track/TrackCmd.h"

namespace PyOpal {
    template <>
    std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<TrackCmd>::attributes = {
        {"LINE", "line", "", PyOpalObjectNS::STRING},
        {"BEAM", "beam", "", PyOpalObjectNS::STRING},
        {"DT", "time_steps", "", PyOpalObjectNS::DOUBLE},  // array
        {"DTSCINIT", "dt_space_charge", "", PyOpalObjectNS::DOUBLE},
        {"DTAU", "dtau", "", PyOpalObjectNS::DOUBLE},
        {"T0", "t0", "", PyOpalObjectNS::DOUBLE},
        {"MAXSTEPS", "max_steps", "", PyOpalObjectNS::FLOAT_LIST},
        {"STEPSPERTURN", "steps_per_turn", "", PyOpalObjectNS::DOUBLE},
        {"ZSTART", "z_start", "", PyOpalObjectNS::DOUBLE},
        {"ZSTOP", "z_stop", "", PyOpalObjectNS::FLOAT_LIST},
        {"TIMEINTEGRATOR", "time_integrator", "", PyOpalObjectNS::PREDEFINED_STRING},
        {"MAP_ORDER", "map_order", "", PyOpalObjectNS::DOUBLE},
    };

    namespace PyTrackCmdNS {

        // Can't use the default PyObject execute function because we need to call
        // setIsParseable to false (otherwise OPAL will try to parse it as an OPAL file)
        void executeWrapper(PyOpalObjectNS::PyOpalObject<TrackCmd>& cmd) {
            std::shared_ptr<TrackCmd> objectPtr = cmd.getOpalShared();
            objectPtr->setIsParseable(false);
            objectPtr->execute();
        }

        BOOST_PYTHON_MODULE(track) {
            ExceptionTranslation::registerExceptions();
            PyOpal::Globals::Initialise();
            PyOpalObjectNS::PyOpalObject<TrackCmd> trackCmd;
            auto trackClass = trackCmd.make_class("Track");
            trackClass.def("execute", &executeWrapper);
        }

    }  // namespace PyTrackCmdNS
}  // namespace PyOpal
