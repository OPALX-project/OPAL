#ifndef __SPLIT_STRATEGY__
#define __SPLIT_STRATEGY__

#include "mpi.h"

#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"


enum commGroupColorings_t {Internal, ExternalToOpt, ExternalToWorker};

/// Defines an interface for splitter strategy implementations
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
            throw OptPilotException("SplitStrategy::SplitStrategy",
                                    "We need 3 or more cores to split!");
    }


    virtual ~SplitStrategy()
    {}


    /**
     *  Forcing concrete implementation to split and assign poller, optimizer
     *  and worker nodes.
     */
    virtual void split() = 0;

    /// Get number of worker groups
    virtual int getNrWorkerGroups() const = 0;

    MPI_Comm getComm()  const { return comm_; }

    int getRank()       const { return rank_; }
    int getGlobalRank() const { return global_rank_; }
    int getNP()         const { return num_procs_; }

    Role_t getRole()    const { return role_; }
    int getLeader()     const { return leader_; }
    int getPoller()     const { return poller_; }


private:

    /// communicator we are splitting
    MPI_Comm comm_;


protected:

    int rank_;
    int global_rank_;
    int num_procs_;
    int group_id_;

    CmdArguments_t cmd_args_;

    Role_t role_;

    /// defines comm splitting
    std::vector<unsigned int> colorings_;

    /// every core specifies a leader (master is its own leader)
    int leader_;

    /// every core can specifies a master
    int poller_;
};

#endif
