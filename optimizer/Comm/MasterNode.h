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

#include <cmath>
#include <set>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include "mpi.h"

#include "Comm/SolutionStateIO.h"


//XXX: SolutionState_t must be serializable, i.e. its value_type must provide
//     writeState(std::ostream&) const and readState(std::istream&)!
//     (call SerializableSolutionState_t?)
template <
      class SolutionState_t
    , class NeighborStrategy_t
>
class MasterNode : public NeighborStrategy_t {

    /// header prefixed to the serialized data in the (single) shared window:
    /// revision together with the actual number of valid payload bytes
    /// (<= buf_size_upper_bound_). Keeping both in the same window as the
    /// payload and writing/reading them within a single fence epoch avoids
    /// tearing between "new revision" and "new data" that two independent
    /// windows could not guarantee.
    struct Meta_t {
        std::size_t revision = 0;
        std::size_t size = 0;
    };

public:
    MasterNode(MPI_Comm master_comm, std::size_t buf_size_upper_bound, std::size_t dim,
               int island_id)
            : buf_size_upper_bound_(buf_size_upper_bound)
            , win_bytes_(sizeof(Meta_t) + buf_size_upper_bound)
            , master_comm_(master_comm)
            , revision_(0) {

        int tmp = 0;
        MPI_Comm_rank(master_comm, &tmp);
        myID_ = static_cast<std::size_t>(tmp);

        MPI_Comm_size(master_comm, &tmp);
        numMasters_ = static_cast<std::size_t>(tmp);
        revision_state_.resize(numMasters_, 0);

        // better to use MPI-2 memory allocation methods
        MPI_Alloc_mem(win_bytes_, MPI_INFO_NULL, &serialized_best_values_);

        // zero-initialize the header so an unwritten window reads as revision 0
        Meta_t empty;
        memcpy(serialized_best_values_, &empty, sizeof(Meta_t));

        // expose our shared memory holding header + best values as one window
        MPI_Win_create(serialized_best_values_, win_bytes_,
                       sizeof(char), MPI_INFO_NULL, master_comm, &win_);

        // execute neighbor strategy to learn which neighbors have to be
        // updated with our solution state (and we collect from)
        collectFrom_ = this->execute(numMasters_, dim, myID_, island_id);
    }

    ~MasterNode() {
        MPI_Win_free(&win_);
        MPI_Free_mem(serialized_best_values_);
    }

    /// store my best values
    void store(char *local_state, std::size_t buffer_size) {

        if (buffer_size > buf_size_upper_bound_) {
            std::cerr << "MasterNode::store(): serialized state (" << buffer_size
                      << " bytes) exceeds window capacity (" << buf_size_upper_bound_
                      << " bytes), dropping this update" << std::endl;
            return;
        }

        Meta_t meta;
        meta.revision = ++revision_;
        meta.size = buffer_size;

        // header and payload are written within the same epoch so a reader
        // can never observe a new revision paired with stale/partial data
        MPI_Win_fence(MPI_MODE_NOPUT, win_);
        memcpy(serialized_best_values_, &meta, sizeof(Meta_t));
        memcpy(serialized_best_values_ + sizeof(Meta_t), local_state, buffer_size);
        MPI_Win_fence(MPI_MODE_NOPUT, win_);
    }

    /// collect all best values from all other masters
    void collect(std::ostringstream &states) {

        char *buffer;
        MPI_Alloc_mem(win_bytes_, MPI_INFO_NULL, &buffer);
        SolutionState_t tmp_states;

        for (std::size_t i=0; i < numMasters_; i++) {
            // ignore all except for selected master PIDs
            if (i == myID_) continue;
            if (collectFrom_.count(i) == 0) continue;

            // cheap pre-check: only the header, to avoid a full fetch when unchanged
            Meta_t hint;
            MPI_Get(&hint, sizeof(Meta_t), MPI_BYTE, i, 0, sizeof(Meta_t), MPI_BYTE, win_);
            MPI_Win_fence(0, win_);

            if(hint.revision <= revision_state_[i]) continue;

            // fetch header + payload together so they can never be torn apart
            MPI_Get(buffer, win_bytes_, MPI_BYTE, i, 0, win_bytes_, MPI_BYTE, win_);
            MPI_Win_fence(0, win_);

            Meta_t meta;
            memcpy(&meta, buffer, sizeof(Meta_t));

            // re-check with the value that is guaranteed consistent with the payload
            if (meta.revision <= revision_state_[i]) continue;
            revision_state_[i] = meta.revision;

            // build the stream from an explicit length, buffer is not null-terminated
            std::istringstream is(std::string(buffer + sizeof(Meta_t), meta.size));

            //XXX: ugly that we have to know the SolutionState_t here
            SolutionState_t state;
            readSolutionState(is, state);
            tmp_states.insert(tmp_states.end(), state.begin(), state.end());
        }

        writeSolutionState(states, tmp_states);

        MPI_Free_mem(buffer);
    }

private:
    /// pointer to MPI window holding header (revision + size) and best solution state
    char *serialized_best_values_;

    /// upper bound on the serialized payload (excludes the header)
    std::size_t buf_size_upper_bound_;
    /// total window size in bytes: sizeof(Meta_t) + buf_size_upper_bound_
    std::size_t win_bytes_;
    std::size_t numMasters_;
    MPI_Comm master_comm_;

    /// single window for header + solution state; all ranks must call
    /// MPI_Win_fence on it the same number of times, in the same order, which
    /// requires collectFrom_ to have equal cardinality on every rank
    MPI_Win win_;

    std::size_t myID_;

    /// neighbors we collect solution states from
    std::set<std::size_t> collectFrom_;
    /// my solution state revision number
    std::size_t revision_;
    /// revision numbers of my neighbors
    std::vector<std::size_t> revision_state_;
};

#endif
