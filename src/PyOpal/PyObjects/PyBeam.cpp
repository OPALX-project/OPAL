//
// Python API for Beam
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

#include "Structure/Beam.h"

namespace PyOpal {
namespace PyBeamNS {

// DOUBLE, STRING, BOOL, INT
template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<Beam>::attributes = {
    {"PARTICLE", "particle", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"MASS", "mass", "", PyOpalObjectNS::DOUBLE},
    {"CHARGE", "charge", "", PyOpalObjectNS::DOUBLE},
    {"ENERGY", "energy", "", PyOpalObjectNS::DOUBLE},
    {"PC", "momentum", "", PyOpalObjectNS::DOUBLE},
    {"GAMMA", "gamma", "", PyOpalObjectNS::DOUBLE},
    {"BCURRENT", "beam_current", "", PyOpalObjectNS::DOUBLE},
    {"BFREQ", "beam_frequency", "", PyOpalObjectNS::DOUBLE},
    {"NPART", "number_of_particles", "", PyOpalObjectNS::DOUBLE},
    {"MOMENTUMTOLERANCE", "momentum_tolerance", "", PyOpalObjectNS::DOUBLE},
};

BOOST_PYTHON_MODULE(beam) {
    ExceptionTranslation::registerExceptions();
    PyOpal::Globals::Initialise();
    PyOpalObjectNS::PyOpalObject<Beam> aBeam;
    auto beamClass = aBeam.make_class("Beam");
    aBeam.addRegister(beamClass);
}

} // PyBeamNS
} // PyOpal

