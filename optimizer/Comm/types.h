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
