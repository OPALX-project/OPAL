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


    virtual void split() {

        parseArguments();

        size_t     group_size  = num_procs_ / num_masters_;
                   group_id_   = rank_ / group_size;
        Comm::id_t group_start = group_id_ * group_size;

        // Fix Pilot to core start_group + 1
        poller_ = group_start;

        // Master and Optimizer fixed to first two cores of group
        if(rank_ % group_size == 0) {

            role_ = POLLER;
            leader_ = 0;

        } else if(rank_ % group_size == 1) {

            role_ = OPTIMIZER;
            leader_ = 1;

        } else {

            role_ = WORKER;
            Comm::localId_t worker_group = ((rank_ % group_size) - 2) /
                                           num_coworkers_worker_;

            leader_ = group_start + 2 + worker_group * num_coworkers_worker_;
            leader_ = leader_ % group_size;
        }

        // define coloring for splitting starting with INTERNAL comm
        colorings_.push_back(group_start + leader_);

        // .. and optimizer -- poller leaders
        if(role_ == WORKER ||
           rank_ % group_size != static_cast<size_t>(leader_))
            colorings_.push_back(MPI_UNDEFINED);
        else
            colorings_.push_back(group_id_);

        // .. and worker -- poller leaders
        if(role_ == OPTIMIZER ||
           rank_ % group_size != static_cast<size_t>(leader_))
            colorings_.push_back(MPI_UNDEFINED);
        else
            colorings_.push_back(group_id_);

        // .. and finally the "world" communicator
        if(role_ == WORKER)
            colorings_.push_back(0);
        else if(role_ == OPTIMIZER)
            colorings_.push_back(1);
        else
            colorings_.push_back(2);

        //FIXME:
        if(role_ == POLLER)
            leader_ = 1;
    }

    virtual int getNrWorkerGroups() const {
        return (num_procs_ - 2 * num_masters_) / num_coworkers_worker_;
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

        if(num_masters_ == 0 || num_coworkers_worker_ == 0) {
            std::cout << "\033[01;31m" << "Need at least"
                      << " 1 master and 1 coworker to run.. Aborting." << "\e[0m" << std::endl;
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
