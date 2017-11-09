#ifndef __SPLIT_STRATEGY__
#define __SPLIT_STRATEGY__

#include "mpi.h"

#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"


enum commGroupColorings_t {Internal, ExternalToOpt, ExternalToWorker};

/**
 * \brief Defines an interface for splitter strategy implementations.
 * \see CommSplitter
 * 
 * A split strategy is a coloring from MPI ranks to colors (ints) to form distinct 
 * communicator groups (MPI_Comms):
 * 
 *   - Internal: all ranks beloning to the same subgroup, i.e. all workers running
 *               one simulation concurrently (as i.e. OPAL needs a "comm group" and 
 *               not a list of ranks).
 *   - ExternalToOpt: the leader ranks to exchange information about the optimization 
 *                    part.
 *   - ExternalToWorker: the leader ranks to start new simulations, get results, ..
 *   - Broadcast: all ranks with the same role to broadcast commands, i.e. STOP.
 * 
 *  Every implementation needs to populate the colorings_ vector, which is turn used 
 *  by the CommSplitter.
 */
class SplitStrategy {

public:

    SplitStrategy(CmdArguments_t args,
                  MPI_Comm comm = MPI_COMM_WORLD)
        : comm_(comm)
        , cmd_args_(args)
        , role_(UNASSIGNED)
    {

        MPI_Comm_rank(comm, &rank_);
        MPI_Comm_rank(MPI_COMM_WORLD, &global_rank_);
        MPI_Comm_size(comm, &num_procs_);

        group_id_ = 0;

        if(num_procs_ < 3)
            throw OptPilotException("SplitStrategy::SplitStartegy",
                                    "We need 3 or more cores to split!");
    }


    virtual ~SplitStrategy()
    {}


    /**
     *  Forcing concrete implementation to split and assign poller, optimizer
     *  and worker nodes.
     */
    virtual void split() = 0;

    MPI_Comm getComm()  const { return comm_; }

    int getRank()           const { return rank_; }
    int getGlobalRank()     const { return global_rank_; }
    int getNP()             const { return num_procs_; }
    int getGroupId()        const { return group_id_; }
    int getWorkerGroupId()  const { return worker_group_id_; }

    Role_t getRole()    const { return role_; }
    int getLeader()     const { return leader_; }
    int getPoller()     const { return poller_; }

    std::vector<int> getWorkers()     const { return workers_; }
    std::vector<int> getOptimizers()  const { return optimizers_; }
    std::vector<int> getCoworkers()   const { return coworkers_; }


private:

    /// communicator we are splitting
    MPI_Comm comm_;


protected:

    int rank_;
    int global_rank_;
    int num_procs_;
    int group_id_;
    int worker_group_id_;

    CmdArguments_t cmd_args_;

    Role_t role_;

    /// defines comm splitting
    std::vector<unsigned int> colorings_;

    /// every core specifies a leader (master is its own leader)
    int leader_;

    /// the leaders local pid of worker groups
    int leader_local_pid_;

    /// the pilot local pid of group
    int master_local_pid_;

    /// the master running the pilot
    int poller_;

    /// used in master <-> workers communicator
    std::vector<int> workers_;

    /// used in aster <-> optimizers communicator
    std::vector<int> optimizers_;

    /// every role has one or more pids to solve the task at hand
    // a worker leader has more worker-cores to solve the forward problem
    // a optimizer leader has more opt-cores to solve the opt problem
    // a master has other masters working on the same opt problem
    std::vector<int> coworkers_;
};

#endif
