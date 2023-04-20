//
// Python API for the FieldSolver
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

#include "Structure/FieldSolver.h"

#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

namespace PyOpal {
namespace PyFieldSolverNS {

// DOUBLE, STRING, BOOL, INT
template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<FieldSolver>::attributes = {
    {"FSTYPE", "type", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"MX", "mesh_size_x", "", PyOpalObjectNS::DOUBLE},
    {"MY", "mesh_size_y", "", PyOpalObjectNS::DOUBLE},
    {"MT", "mesh_size_t", "", PyOpalObjectNS::DOUBLE},
    {"PARFFTX", "parallelise_fft_x", "", PyOpalObjectNS::BOOL},
    {"PARFFTY", "parallelise_fft_y", "", PyOpalObjectNS::BOOL},
    {"PARFFTT", "parallelise_fft_t", "", PyOpalObjectNS::BOOL},
    {"BCFFTX", "fft_boundary_x", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"BCFFTY", "fft_boundary_y", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"BCFFTZ", "fft_boundary_z", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"GREENSF", "greens_function", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"BBOXINCR", "bounding_box_increase", "", PyOpalObjectNS::DOUBLE},
    {"GEOMETRY", "geometry", "", PyOpalObjectNS::UPPER_CASE_STRING},
    {"ITSOLVER", "iterative_solver", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"INTERPL", "interpolation", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"TOL", "tolerance", "", PyOpalObjectNS::DOUBLE},
    {"MAXITERS", "max_iterations", "", PyOpalObjectNS::DOUBLE},
    {"PRECMODE", "preconditioner_mode", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"RC", "cutoff_radius", "", PyOpalObjectNS::DOUBLE},
    {"ALPHA", "alpha", "", PyOpalObjectNS::DOUBLE},
};

void registerFieldSolver(PyOpalObjectNS::PyOpalObject<FieldSolver>& pyfs) {
    Object* obj = &(*pyfs.getOpalShared());
    FieldSolver* fs = dynamic_cast<FieldSolver*>(obj);
    if (fs == nullptr) {
        throw OpalException(
            "PyOpal::PyFieldSolverNS::registerFieldSolver",
            "Internal error - field solver not recognised during register()"
        );
    }
    fs->execute();
    OpalData::getInstance()->define(obj);
}

BOOST_PYTHON_MODULE(field_solver) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<FieldSolver> fs;
    auto fsClass = fs.make_class("FieldSolver");
    fsClass.def("register", &registerFieldSolver);
}

} // PyFieldSolverNS
} // PyOpal

