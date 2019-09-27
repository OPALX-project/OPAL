#ifndef __WORKER_H__
#define __WORKER_H__

#include <iostream>

#include "boost/smart_ptr.hpp"

#include "Pilot/Poller.h"

#include "Comm/types.h"
#include "Util/Types.h"
#include "Util/MPIHelper.h"
#include "Util/CmdArguments.h"

/**
 *  \class Worker
 *  \brief A worker MPI entity consists of a processor group that runs a
 *  simulation of type Sim_t. The main loop in run() accepts new jobs from the
 *  master process runs the simulation and reports back the results.
 *
 *  @see Pilot
 *  @see Poller
 *  @see MPIHelper.h
 *
 *  @tparam Sim_T type of simulation to run
 */
template <class Sim_t>
class Worker : protected Poller {

public:

    Worker(Expressions::Named_t constraints,
           std::string simName,
           Comm::Bundle_t comms,
           CmdArguments_t args)
        : Poller(comms.worker)
        , cmd_args_(args)
    {
        constraints_     = constraints;
        simulation_name_ = simName;
        pilot_rank_      = comms.master_local_pid;
        is_idle_         = true;
        coworker_comm_   = comms.coworkers;

        leader_pid_      = 0;
        MPI_Comm_size(coworker_comm_, &num_coworkers_);
    }

    Worker(Expressions::Named_t objectives,
           Expressions::Named_t constraints,
           std::string simName,
           Comm::Bundle_t comms,
           CmdArguments_t args,
           bool isOptimizer = true,
           std::vector<std::string> filesToKeep = {})
        : Poller(comms.worker)
        , cmd_args_(args)
        , filesToKeep_m(filesToKeep)
    {
        objectives_      = objectives;
        constraints_     = constraints;
        simulation_name_ = simName;
        pilot_rank_      = comms.master_local_pid;
        is_idle_         = true;
        coworker_comm_   = comms.coworkers;

        leader_pid_      = 0;
        MPI_Comm_size(coworker_comm_, &num_coworkers_);

        if (!isOptimizer) return;
        int my_local_pid = 0;
        MPI_Comm_rank(coworker_comm_, &my_local_pid);

        // distinction between leader and coworkers
        if(my_local_pid == leader_pid_)
            run();
        else
            runCoWorker();
    }

    ~Worker()
    {}


protected:
    typedef boost::scoped_ptr<Sim_t> SimPtr_t;

    bool is_idle_;
    MPI_Comm coworker_comm_;

    Expressions::Named_t objectives_;
    Expressions::Named_t constraints_;


    /// coworkers simply wait on a job broadcast from the leader and then
    /// start a simulation..
    void runCoWorker() {

        MPI_Request stop_req;
        size_t stop_value = 0;

        MPI_Irecv(&stop_value, 1, MPI_UNSIGNED_LONG, leader_pid_,
                  MPI_ANY_TAG, coworker_comm_, &stop_req);
        is_running_ = true;

        while(is_running_) {

            //FIXME: bcast blocks after our leader stopped working
            // Either we create a new class implementing a coworker in the
            // same manner as the worker (poll loop). Anyway there is no way
            // around removing the Bcast and adding another tag in the poll
            // loop above in order to be able to exit cleanly.
            if(stop_req != MPI_REQUEST_NULL) {
                MPI_Status status;
                int flag = 0;
                MPI_Test(&stop_req, &flag, &status);

                if(flag) {

                    if(status.MPI_TAG == MPI_COWORKER_NEW_JOB_TAG) {
                        Param_t params;
                        MPI_Bcast_params(params, leader_pid_, coworker_comm_);

                        try {
                            SimPtr_t sim(new Sim_t(objectives_, constraints_,
                                    params, simulation_name_, coworker_comm_,
                                    cmd_args_));

                            sim->run();
                        } catch(OptPilotException &ex) {
                            std::cout << "Exception while running simulation: "
                                      << ex.what() << std::endl;
                        }
                        MPI_Irecv(&stop_value, 1, MPI_UNSIGNED_LONG, leader_pid_,
                                  MPI_ANY_TAG, coworker_comm_, &stop_req);
                    }

                    if(status.MPI_TAG == MPI_STOP_TAG) {
                        is_running_ = false;
                        break;
                    }
                }
            }
        }
    }


protected:

    int leader_pid_;
    int num_coworkers_;
    int pilot_rank_;
    std::string simulation_name_;
    CmdArguments_t cmd_args_;

    /// notify coworkers of incoming broadcast
    void notifyCoWorkers(int tag) {

        for(int i=0; i < num_coworkers_; i++) {
            if(i == leader_pid_) continue;

            size_t dummy = 0;
            MPI_Send(&dummy, 1, MPI_UNSIGNED_LONG, i, tag, coworker_comm_);
        }
    }

    void setupPoll() {
        size_t dummy = 1;
        MPI_Send(&dummy, 1, MPI_UNSIGNED_LONG, pilot_rank_,
                 MPI_WORKER_STATUSUPDATE_TAG, comm_m);
    }

    void prePoll()
    {}

    void postPoll()
    {}

    void onStop() {
        if(num_coworkers_ > 1)
            notifyCoWorkers(MPI_STOP_TAG);
    }

    virtual bool onMessage(MPI_Status status, size_t recv_value) {

        if(status.MPI_TAG == MPI_WORK_JOBID_TAG) {

            is_idle_ = false;
            size_t job_id = recv_value;

            // get new job
            Param_t params;
            MPI_Recv_params(params, (size_t)pilot_rank_, comm_m);

            // and forward to coworkers (if any)
            if(num_coworkers_ > 1) {
                notifyCoWorkers(MPI_COWORKER_NEW_JOB_TAG);
                MPI_Bcast_params(params, leader_pid_, coworker_comm_);
            }

            //XXX we need to know if we want EVAL or DERIVATIVE
            //reqVarContainer_t reqVars;
            //MPI_Recv_reqvars(reqVars, (size_t)pilot_rank_, comm_m);

            reqVarContainer_t requested_results;
            std::string simdir;
            try {
                SimPtr_t sim(new Sim_t(objectives_, constraints_,
                        params, simulation_name_, coworker_comm_, cmd_args_));

                // run simulation in a "blocking" fashion
                sim->run();
                sim->collectResults();

                if ( filesToKeep_m.empty() ) {
                    sim->cleanUp();
                } else {
                    // if empty, we keep all files
                    sim->cleanUp(filesToKeep_m);
                }

                requested_results = sim->getResults();
                simdir = sim->getSimDirectory();

            } catch(OptPilotException &ex) {
                std::cout << "Exception while running simulation: "
                          << ex.what() << std::endl;
            }

            MPI_Send(&job_id, 1, MPI_UNSIGNED_LONG, pilot_rank_,
                     MPI_WORKER_FINISHED_TAG, comm_m);

            MPI_Send(simdir.c_str(), simdir.size(), MPI_CHAR, pilot_rank_,
                     MPI_WORKER_DIRECTORY_TAG, comm_m);

            size_t dummy = 0;
            MPI_Recv(&dummy, 1, MPI_UNSIGNED_LONG, pilot_rank_,
                     MPI_WORKER_FINISHED_ACK_TAG, comm_m, &status);

            MPI_Send_reqvars(requested_results, (size_t)pilot_rank_, comm_m);


            is_idle_ = true;
            return true;

        } else {
            std::stringstream os;
            os << "Unexpected MPI_TAG: " << status.MPI_TAG;
            std::cout << "(Worker) Error: " << os.str() << std::endl;
            throw OptPilotException("Worker::onMessage", os.str());
        }
    }
private:
    const std::vector<std::string> filesToKeep_m;
};

#endif