#ifndef __SAMPLE_WORKER_H__
#define __SAMPLE_WORKER_H__

#include "Pilot/Worker.h"


/**
 *  \class SampleWorker
 *  \brief A worker MPI entity consists of a processor group that runs a
 *  simulation of type Sim_t. The main loop in run() accepts new jobs from the
 *  master process runs the simulation and reports back the results.
 *
 *  @see SamplePilot
 *  @see Worker
 *  @see MPIHelper.h
 *
 *  @tparam Sim_T type of simulation to run
 */
template <class Sim_t>
class SampleWorker : protected Worker<Sim_t> {

public:

    SampleWorker(Expressions::Named_t objectives,
                 Expressions::Named_t constraints,
                 std::string simName,
                 Comm::Bundle_t comms,
                 CmdArguments_t args,
                 const std::vector<std::string> &storeobjstr,
                 const std::vector<std::string> &filesToKeep,
                 const std::map<std::string, std::string> &userVariables)
        : Worker<Sim_t>(objectives, constraints, simName, comms, args, userVariables, false)
        , statVariablesToStore_m(storeobjstr)
        , filesToKeep_m(filesToKeep)
    {

        int my_local_pid = 0;
        MPI_Comm_rank(this->coworker_comm_, &my_local_pid);

        // distinction between leader and coworkers
        if(my_local_pid == this->leader_pid_)
            this->run();
        else
            runSlave();
    }

    ~SampleWorker()
    {}

protected:
    // FIXME Open issue #250 (https://gitlab.psi.ch/OPAL/src/issues/250)
    const std::vector<std::string> statVariablesToStore_m;

    /// notify coworkers of incoming broadcast
    void notifyCoWorkers(size_t job_id, int tag) {

        for(int i=0; i < this->num_coworkers_; i++) {
            if(i == this->leader_pid_) continue;

            // send job id to co workers
            MPI_Send(&job_id, 1, MPI_UNSIGNED_LONG, i, tag, this->coworker_comm_);
        }
    }

    /// coworkers simply wait on a job broadcast from the leader and then
    /// start a simulation..
    void runSlave() {
        /* needs to be executed by derived class otherwise
         * a base class instance is created.
         */

        MPI_Request stop_req;
        size_t job_id = 0;

        MPI_Irecv(&job_id, 1, MPI_UNSIGNED_LONG, this->leader_pid_,
                  MPI_ANY_TAG, this->coworker_comm_, &stop_req);
        this->is_running_ = true;

        while(this->is_running_) {

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
                        MPI_Bcast_params(params, this->leader_pid_, this->coworker_comm_);

                        try {
                            typename Worker<Sim_t>::SimPtr_t sim(
                                new Sim_t(this->objectives_, this->constraints_,
                                          params, this->simulation_name_, this->coworker_comm_,
                                          this->cmd_args_, this->userVariables_));

                            sim->setFilename(job_id);

                            sim->run();
                        } catch(OptPilotException &ex) {
                            std::cout << "Exception while running simulation: "
                                      << ex.what() << std::endl;
                        }
                        MPI_Irecv(&job_id, 1, MPI_UNSIGNED_LONG, this->leader_pid_,
                                  MPI_ANY_TAG, this->coworker_comm_, &stop_req);
                    }

                    if(status.MPI_TAG == MPI_STOP_TAG) {
                        this->is_running_ = false;
                        break;
                    }
                }
            }
        }
    }

    bool onMessage(MPI_Status status, size_t recv_value) override {

        if(status.MPI_TAG == MPI_WORK_JOBID_TAG) {

            this->is_idle_ = false;
            size_t job_id = recv_value;

            // get new job
            Param_t params;
            MPI_Recv_params(params, (size_t)this->pilot_rank_, this->comm_m);

            // and forward to coworkers (if any)
            if(this->num_coworkers_ > 1) {
                notifyCoWorkers(job_id, MPI_COWORKER_NEW_JOB_TAG);
                MPI_Bcast_params(params, this->leader_pid_, this->coworker_comm_);
            }

            reqVarContainer_t requested_results;
            try {
                typename Worker<Sim_t>::SimPtr_t sim(new Sim_t(this->objectives_,
                                                               this->constraints_,
                                                               params,
                                                               this->simulation_name_,
                                                               this->coworker_comm_,
                                                               this->cmd_args_,
                                                               this->userVariables_));

                sim->setFilename(job_id);

                // run simulation in a "blocking" fashion
                sim->run();

                // this is requests the columns from the stat file and stores them
                // in a map with the column names as key and the columns as values; for #250
                //
                // std::map<std::string,
                //          std::vector<double> > data = sim->getData(statVariablesToStore_m);

                sim->collectResults();
                requested_results = sim->getResults();

                // base clase of SamplePilot requires at least 1 objective --> dummy objective (SamplePilot, line 64)
                if ( (this->objectives_.size() > 1) && filesToKeep_m.empty() ) {
                    sim->cleanUp();
                } else {
                    // if empty, we keep all files
                    sim->cleanUp(filesToKeep_m);
                }

            } catch(OptPilotException &ex) {
                std::cout << "Exception while running simulation: "
                          << ex.what() << std::endl;
            }

            MPI_Send(&job_id, 1, MPI_UNSIGNED_LONG, this->pilot_rank_,
                     MPI_WORKER_FINISHED_TAG, this->comm_m);

            size_t dummy = 0;
            MPI_Recv(&dummy, 1, MPI_UNSIGNED_LONG, this->pilot_rank_,
                     MPI_WORKER_FINISHED_ACK_TAG, this->comm_m, &status);

            MPI_Send_reqvars(requested_results, (size_t)this->pilot_rank_, this->comm_m);

            this->is_idle_ = true;
            return true;

        } else {
            std::stringstream os;
            os << "Unexpected MPI_TAG: " << status.MPI_TAG;
            std::cout << "(Worker) Error: " << os.str() << std::endl;
            throw OptPilotException("SampleWorker::onMessage", os.str());
        }
    }

private:
    const std::vector<std::string> filesToKeep_m;
};

#endif