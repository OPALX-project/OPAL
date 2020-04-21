//
// Types in namespace Comm
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
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
#ifndef __COMM_TYPES__
#define __COMM_TYPES__

#include <vector>

#include "mpi.h"

namespace Comm {

    typedef size_t id_t;
    typedef size_t localId_t;

    /// bundles all communicators for a specific role/pid
    struct Bundle_t {
        int island_id;
        int leader_pid;
        int master_pid;
        int master_local_pid;
        MPI_Comm worker;
        MPI_Comm opt;
        MPI_Comm coworkers;
        MPI_Comm world;
    };
}

#endif
