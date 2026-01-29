//
// Class MasterNode
//   Implements a node in the network of all pilots, exposing store and
//   collect operations on a specific set of neighbors.
//
//   Using the neighbor strategy a set of neighbors we collect solution state
//   from (and they collect from us) is defined. Using this set of neighbors the
//   solution states propagate throughout the network. The store and collect
//   operations are implemented using one sided MPI communication methods
//   (simulating shared memory).
//   A revision number is used to prevent receiving previously collected
//   solution states from neighbors.
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
#ifndef __MASTER_NODE__
#define __MASTER_NODE__

#include <set>
#include <cmath>
#include <vector>

#include <string>
#include <sstream>

#include "mpi.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>


//XXX: SolutionState_t must be serializable! (call
//     SerializableSolutionState_t?)
template <
      class SolutionState_t
    , class NeighborStrategy_t
>
class MasterNode : public NeighborStrategy_t {


public:

    MasterNode(MPI_Comm master_comm, size_t buf_size_upper_bound, size_t dim,
               int island_id)
            : buf_size_upper_bound_(buf_size_upper_bound)
            , master_comm_(master_comm) {

        int tmp = 0;
        MPI_Comm_rank(master_comm, &tmp);
        myID_ = static_cast<size_t>(tmp);

        MPI_Comm_size(master_comm, &tmp);
        numMasters_ = static_cast<size_t>(tmp);
        revision_state_.resize(numMasters_, 0);

        // better to use MPI-2 memory allocation methods
        MPI_Alloc_mem(sizeof(char) * buf_size_upper_bound,
                      MPI_INFO_NULL, &serialized_best_values_);

        // expose our shared memory holding our best values
        MPI_Win_create(serialized_best_values_, buf_size_upper_bound,
                       sizeof(char), MPI_INFO_NULL, master_comm, &win_);

        MPI_Win_create(&revision_, 1, sizeof(size_t), MPI_INFO_NULL,
                       master_comm, &win_rev_);


        // execute neighbor strategy to learn which neighbors have to be
        // updated with our solution state (and we collect from)
        collectFrom_ = this->execute(numMasters_, dim, myID_, island_id);
    }

    ~MasterNode() {
        MPI_Win_free(&win_);
        MPI_Free_mem(serialized_best_values_);
        MPI_Win_free(&win_rev_);
    }


    /// store my best values
    void store(char *local_state, size_t buffer_size) {

        revision_++;
        MPI_Win_fence(MPI_MODE_NOPUT, win_rev_);

        size_t buf_size = buffer_size;
        if(buf_size > buf_size_upper_bound_)
            std::cout << "windows too small: " << buffer_size << " / "
                      << buf_size_upper_bound_ << std::endl;

        memcpy(serialized_best_values_, local_state, buf_size);
        MPI_Win_fence(MPI_MODE_NOPUT, win_);
        //MPI_Win_fence((MPI_MODE_NOPUT | MPI_MODE_NOPRECEDE), win_);
    }


    /// collect all best values from all other masters
    void collect(std::ostringstream &states) {

        char *buffer;
        MPI_Alloc_mem(buf_size_upper_bound_, MPI_INFO_NULL, &buffer);
        SolutionState_t tmp_states;


        for(size_t i=0; i < numMasters_; i++) {
            // ignore all except for selected master PIDs
            if(i == myID_) continue;
            if(collectFrom_.count(i) == 0) continue;

            // only continue if new values are available on master i
            size_t revision = 0;
            MPI_Get(&revision, 1, MPI_UNSIGNED_LONG, i, 0, 1, MPI_UNSIGNED_LONG, win_rev_);
            MPI_Win_fence(0, win_rev_);

            if(revision <= revision_state_[i]) continue;
            revision_state_[i] = revision;

            MPI_Get(buffer, buf_size_upper_bound_, MPI_CHAR, i, 0,
                    buf_size_upper_bound_, MPI_CHAR, win_);
            MPI_Win_fence(0, win_);

            std::istringstream is(buffer);
            boost::archive::text_iarchive ia(is);

            //XXX: ugly that we have to know the SolutionState_t here
            SolutionState_t state;
            ia >> state;
            tmp_states.insert(tmp_states.end(), state.begin(), state.end());
        }

        boost::archive::text_oarchive oa(states);
        oa << tmp_states;

        MPI_Free_mem(buffer);
    }


private:

    /// pointer to MPI window holding current best solution state
    char *serialized_best_values_;

    /// and upper bound on the allocated memory in the MPI window
    size_t buf_size_upper_bound_;
    size_t numMasters_;
    MPI_Comm master_comm_;

    // windows for storing revision and solution state
    MPI_Win win_;
    MPI_Win win_rev_;

    size_t myID_;

    /// neighbors we collect solution states from
    std::set<size_t> collectFrom_;
    /// my solution state revision number
    size_t revision_;
    /// revision numbers of my neighbors
    std::vector<size_t> revision_state_;
};

#endif
