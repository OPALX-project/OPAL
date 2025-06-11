//
// Python API for Option
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

#include "BasicActions/Option.h"

namespace PyOpal {
    // DOUBLE, STRING, BOOL, INT
    template <>
    std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<Option>::attributes = {
        {"ECHO", "echo", "", PyOpalObjectNS::BOOL},
        {"INFO", "info", "", PyOpalObjectNS::BOOL},
        {"TRACE", "trace", "", PyOpalObjectNS::BOOL},
        {"WARN", "warn", "", PyOpalObjectNS::BOOL},
        {"TELL", "tell", "", PyOpalObjectNS::BOOL},
        {"SEED", "seed", "", PyOpalObjectNS::DOUBLE},
        {"PSDUMPFREQ", "ps_dump_frequency", "", PyOpalObjectNS::DOUBLE},
        {"STATDUMPFREQ", "stat_dump_frequency", "", PyOpalObjectNS::DOUBLE},
        {"SCSOLVEFREQ", "sc_solve_frequency", "", PyOpalObjectNS::DOUBLE},
        {"SPTDUMPFREQ", "spt_dump_frequency", "", PyOpalObjectNS::DOUBLE},
        {"MTSSUBSTEPS", "mts_substeps", "", PyOpalObjectNS::DOUBLE},
        {"REMOTEPARTDEL", "remote_particle_delete", "", PyOpalObjectNS::DOUBLE},
        {"PSDUMPFRAME", "ps_dump_frame", "", PyOpalObjectNS::PREDEFINED_STRING},
        {"REPARTFREQ", "repartition_frequency", "", PyOpalObjectNS::DOUBLE},
        {"MINBINEMITTED", "min_bin_emitted", "", PyOpalObjectNS::DOUBLE},
        {"MINSTEPFORREBIN", "min_step_for_rebin", "", PyOpalObjectNS::DOUBLE},
        {"REBINFREQ", "rebin_frequency", "", PyOpalObjectNS::DOUBLE},
        {"RHODUMP", "rho_dump", "", PyOpalObjectNS::BOOL},
        {"EBDUMP", "eb_dump", "", PyOpalObjectNS::BOOL},
        {"CSRDUMP", "csr_dump", "", PyOpalObjectNS::BOOL},
        {"AUTOPHASE", "autophase", "", PyOpalObjectNS::DOUBLE},
        {"CZERO", "czero", "", PyOpalObjectNS::BOOL},
        {"RNGTYPE", "rng_type", "", PyOpalObjectNS::PREDEFINED_STRING},
        {"CLOTUNEONLY", "clo_tune_only", "", PyOpalObjectNS::BOOL},
        {"NUMBLOCKS", "num_blocks", "", PyOpalObjectNS::DOUBLE},
        {"RECYCLEBLOCKS", "recycle_blocks", "", PyOpalObjectNS::DOUBLE},
        {"NLHS", "nlhs", "", PyOpalObjectNS::DOUBLE},
        {"ENABLEHDF5", "enable_hdf5", "", PyOpalObjectNS::BOOL},
        {"ENABLEVTK", "enable_vtk", "", PyOpalObjectNS::BOOL},
        {"ASCIIDUMP", "ascii_dump", "", PyOpalObjectNS::BOOL},
        {"BOUNDPDESTROYFQ", "bound_p_destroy", "", PyOpalObjectNS::DOUBLE},
        {"BEAMHALOBOUNDARY", "beam_halo_boundary", "", PyOpalObjectNS::DOUBLE},
        {"IDEALIZED", "idealized", "", PyOpalObjectNS::BOOL},
        {"LOGBENDTRAJECTORY", "log_bend_trajectory", "", PyOpalObjectNS::BOOL},
        {"VERSION", "version", "", PyOpalObjectNS::DOUBLE}};

    namespace PyOptionNS {

        BOOST_PYTHON_MODULE(option) {
            ExceptionTranslation::registerExceptions();
            PyOpal::Globals::Initialise();
            PyOpalObjectNS::PyOpalObject<Option> anOption;
            auto optionClass = anOption.make_class("Option");
            anOption.addExecute(optionClass);
        }

    }  // namespace PyOptionNS
}  // namespace PyOpal
