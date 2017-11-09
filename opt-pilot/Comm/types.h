#ifndef __COMM_TYPES__
#define __COMM_TYPES__

#include <vector>
#include <map>
#include <set>

#include "mpi.h"

#include "boost/tuple/tuple.hpp"
#include "boost/variant/variant.hpp"

namespace Comm {

    typedef size_t id_t;
    typedef size_t localId_t;

    // defining different groups of processors (consecutive vs. sets).. This
    // is currently not used.
    typedef boost::tuple<size_t, size_t> blockProcessorGroup_t;
    typedef std::set<size_t>             setProcessorGroup_t;

    typedef boost::variant <
            std::vector<blockProcessorGroup_t>
            , std::vector<setProcessorGroup_t>
            >
        processorGroups_t;

    /// bundles all communicators for a specific role/pid
    struct Bundle_t {
        /// unique island identifier (group of pilot, opt, workers)
        int island_id;

        /// global pid of the group (pilot, opt, worker) leader
        int leader_pid;

        /// global pid of the pilot for this group
        int master_pid;

        /// local pid of the pilot in the communicator groups
        int master_local_pid;

        /// communicator of all workers
        MPI_Comm worker;

        /// communicator of all optimizers
        MPI_Comm opt;

        /// communicator of coworkers, i.e. a subset of ranks assigned to one worker
        MPI_Comm coworkers;

        /// world communicator as passed to the CommSplitter
        MPI_Comm world;
    };
}

#endif
