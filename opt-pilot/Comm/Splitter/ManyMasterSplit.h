#ifndef __ONE_MASTER_SPLIT__
#define __ONE_MASTER_SPLIT__

#include "mpi.h"

#include "Comm/types.h"
#include "Comm/Splitter/SplitStrategy.h"

#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"


/**
 *  A very simple splitting strategy where we have a one core optimizer and
 *  pilot (and k of those "islands") and many-core worker groups.
 *  The number of islands and co-workers is retrieved from the passed command
 *  line arguments:
 *    - num-masters
 *    - num-coworkers
 */
template < class TopoDiscoveryStrategy >
class ManyMasterSplit : protected SplitStrategy, public TopoDiscoveryStrategy {

public:

    ManyMasterSplit(CmdArguments_t args, MPI_Comm comm = MPI_COMM_WORLD)
        : SplitStrategy(args, comm)
    {}


    virtual ~ManyMasterSplit()
    {}


    void split() {

        parseArguments();

        size_t     group_size  = num_procs_ / num_masters_;
                   group_id_   = rank_ / group_size;
        Comm::id_t group_start = group_id_ * group_size;

        // fix Pilot to the of the group
        poller_ = group_start;
        master_local_pid_ = 0;

        worker_group_id_ = 0;

        // Pilot/Master and Optimizer fixed to first two cores of group
        if(rank_ % group_size == 0) {

            role_ = POLLER;
            leader_local_pid_ = 0;

        } else if(rank_ % group_size == 1) {

            role_ = OPTIMIZER;
            leader_local_pid_ = 1;

        } else {

            role_ = WORKER;
            worker_group_id_ = ((rank_ - 2) % (group_size - 2)) /
                                           num_coworkers_worker_;

            leader_local_pid_ = 2 + worker_group_id_ * num_coworkers_worker_;

            // worker group id is global rank of leader
            worker_group_id_ = group_start + leader_;
        }

        // define coloring for splitting starting with INTERGROUP comm
        // this is always the "color" of the leader (global rank)
        leader_ = leader_local_pid_ + group_start;
        colorings_.push_back(leader_);

        // now we can create individual comm group between leaders of
        // different groups, i.e. the optimizer leader and pilot leader.
        // all non-participating ranks put MPI_UNDEFINED.

        // for optimizer -- pilot leaders
        if(role_ == WORKER ||
           rank_ % group_size != static_cast<size_t>(local_leader_pid_))
            colorings_.push_back(MPI_UNDEFINED);
        else
            colorings_.push_back(group_id_);

        // for worker -- pilot leaders
        if(role_ == OPTIMIZER ||
           rank_ % group_size != static_cast<size_t>(local_leader_pid_))
            colorings_.push_back(MPI_UNDEFINED);
        else
            colorings_.push_back(worker_group_id_);

        // and finally the "world" communicator for all ranks 
        // with the same role
        if(role_ == WORKER)
            colorings_.push_back(0);
        else if(role_ == OPTIMIZER)
            colorings_.push_back(1);
        else
            colorings_.push_back(2);

        //FIXME: pilot/master is its own leader? unused most likely!
        if(role_ == POLLER) {
            poller_ = MPI_UNDEFINED;
            local_leader_pid_ = MPI_UNDEFINED;
        }
    }

private:

    size_t num_masters_;
    size_t num_coworkers_worker_;

    void parseArguments() {

        num_coworkers_worker_ = 0;
        try {
            num_coworkers_worker_ = cmd_args_->getArg<size_t>("num-coworkers");
        } catch (OptPilotException &e) {
            std::cout << "\033[01;31m" << "Could not find 'num-coworkers' "
                      << "in arguments.. Aborting." << "\e[0m" << std::endl;
            MPI_Abort(getComm(), -111);
        }

        num_masters_ = 1;
        try {
            num_masters_ = cmd_args_->getArg<size_t>("num-masters");
        } catch (OptPilotException &e) {
            std::cout << "\033[01;31m" << "Could not find 'num-masters' "
                      << "in arguments.. Aborting." << "\e[0m" << std::endl;
            MPI_Abort(getComm(), -1111);
        }

        if(static_cast<size_t>(num_procs_) < num_masters_ *
                                             (2 + num_coworkers_worker_)) {
            std::cout << "\033[01;31m" << "Need at least "
                      << (num_coworkers_worker_ + 2) * num_masters_
                      << " cores to run.. Aborting." << "\e[0m" << std::endl;
            MPI_Abort(getComm(), -1111);
        }
    }

};

#endif
