#ifndef __COMM_SPLITTER__
#define __COMM_SPLITTER__

#include "mpi.h"

#include "Comm/types.h"
#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"

//TODO: what is the performance difference between using MPI_COMM_WORLD and
//      p2p communication vs. communicator groups??
/**
 *  \brief Role assignment according to strategy (that might use hardware
 *         network information).
 *
 *  The CommSplitter splits the passed (usually MPI_COMM_WORLD) communicator
 *  into several comm groups using the colors provided by the splitting
 *  strategy.
 *  After construction each processor has an assigned role (optimizer, worker
 *  or pilot) and a set of communicators to send and receive tasks.
 *  The expected colors have the following meaning:
 *
 *    - color[0] is shared by all my co-workers
 *    - color[1] is shared by the optimizer leader and the pilot
 *    - color[2] is shared by the worker leader and the pilot
 *    - color[3] is shared by all processors with the same role (used for
 *      broadcasts)
 *
 */
template< class Strategy_t >
class CommSplitter : public Strategy_t {

public:

    CommSplitter(CmdArguments_t args, MPI_Comm comm = MPI_COMM_WORLD)
        : Strategy_t(args, comm)
        , world_comm_(comm)
    {
        MPI_Comm_group(world_comm_, &world_group_);
        MPI_Comm_rank(world_comm_,  &global_rank_);

        // the splitter strategy computes the colorings
        Strategy_t::split();

        my_opt_comm_      = MPI_COMM_NULL;
        my_worker_comm_   = MPI_COMM_NULL;
        my_coworker_comm_ = MPI_COMM_NULL;


        MPI_Comm_split(world_comm_, Strategy_t::colorings_[0],
                       global_rank_, &my_coworker_comm_ );

        MPI_Comm_split(world_comm_, Strategy_t::colorings_[1],
                       global_rank_, &my_opt_comm_ );

        MPI_Comm_split(world_comm_, Strategy_t::colorings_[2],
                       global_rank_, &my_worker_comm_ );

        MPI_Comm_split(world_comm_, Strategy_t::colorings_[3],
                       global_rank_, &my_comm_world_ );

        // just a precaution to make sure everybody is ready
        MPI_Barrier(world_comm_);
    }

    virtual ~CommSplitter()
    {
        MPI_Group_free(&world_group_);
        if(world_comm_ != MPI_COMM_WORLD)
            MPI_Comm_free(&world_comm_);
    }

    bool isOptimizer() const { return Strategy_t::role_ == OPTIMIZER; }
    bool isWorker()    const { return Strategy_t::role_ == WORKER; }
    bool isPilot()     const { return Strategy_t::role_ == POLLER; }

    MPI_Comm mpiComm() const { return world_comm_; }
    int globalRank()   const { return global_rank_; }
    int pilotRank()    const { return Strategy_t::poller_; }
    int getLeader()    const { return Strategy_t::leader_; }

    int getNrWorkerGroups() const { return Strategy_t::getNrWorkerGroups();}

    /// construct comm bundle and return
    Comm::Bundle_t getBundle() const {

        Comm::Bundle_t bundle;

        bundle.island_id        = Strategy_t::group_id_;

        bundle.leader_pid       = Strategy_t::leader_;
        bundle.master_pid       = Strategy_t::poller_;
        //FIXME: is it always 0?
        bundle.master_local_pid = 0; //poller_local_pid_;

        bundle.worker           = my_worker_comm_;
        bundle.opt              = my_opt_comm_;
        bundle.coworkers        = my_coworker_comm_;
        bundle.world            = my_comm_world_;

        return bundle;
    }

private:

    MPI_Comm  world_comm_;
    MPI_Group world_group_;

    /// global MPI PID
    int global_rank_;

    // define various communicators for communication between components
    MPI_Comm my_worker_comm_;
    MPI_Comm my_opt_comm_;
    MPI_Comm my_coworker_comm_;
    MPI_Comm my_comm_world_;

    /// local (wrt. the communicator group) rank of the master/pilot process
    int poller_local_pid_;
};

#endif
