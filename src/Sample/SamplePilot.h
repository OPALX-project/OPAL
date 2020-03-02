#ifndef __SAMPLE_PILOT_H__
#define __SAMPLE_PILOT_H__

#include "Pilot/Pilot.h"
#include "Sample/SampleWorker.h"
#include "Expression/Parser/function.hpp"

/**
 *  \class SamplePilot
 *  \brief The sample Pilot (Master): Coordinates requests by sampler
 *  to workers.
 *
 *  Every worker thread notifies the master here if idle or not. When
 *  available the master dispatches one of the pending simulations to the
 *  worker who will run the specified simulation and report results back to
 *  the master.
 *
 *  @see SampleWorker
 *  @see Sampler
 *
 *  @tparam Opt_t type of the sampler
 *  @tparam Sim_t type of the simulation
 *  @tparam SolPropagationGraph_t strategy to distribute solution between
 *          master islands
 *  @tparam Comm_t comm splitter strategy
 */
template <
      class Opt_t
    , class Sim_t
    , class SolPropagationGraph_t
    , class Comm_t
    >
class SamplePilot : protected Pilot<Opt_t,
                                    Sim_t,
                                    SolPropagationGraph_t,
                                    Comm_t>
{

public:

    SamplePilot(CmdArguments_t args, boost::shared_ptr<Comm_t> comm,
                functionDictionary_t known_expr_funcs,
                const DVarContainer_t &dvar,
                const Expressions::Named_t &obj,
                const std::map< std::string,
                                std::shared_ptr<SamplingMethod>
                              >& sampleMethods,
                const std::vector<std::string> &storeobjstr,
                const std::vector<std::string> &filesToKeep,
                const std::map<std::string, std::string> &userVariables)
        : Pilot<Opt_t,
                Sim_t,
                SolPropagationGraph_t,
                Comm_t>(args,
                        comm,
                        known_expr_funcs,
                        dvar,
                        obj,
                        Expressions::Named_t(),
                        {},
                        false,
                        {})
        , sampleMethods_m(sampleMethods)
    {
        if (obj.size() == 0) {
        // create a dummy objective, base class requires at least 1 objective
            this->objectives_ = {
                {"dummy", new Expressions::Expr_t("dummy")}
            };
        }

        this->setup(known_expr_funcs, storeobjstr, filesToKeep, userVariables);
    }

    virtual ~SamplePilot()
    {}


protected:

    /// keep track of requests and running jobs
    typedef std::map<size_t, Param_t > Jobs_t;
    typedef Jobs_t::iterator JobIter_t;
    Jobs_t  running_job_list_;
    Jobs_t  request_queue_;


    virtual
    void setup(functionDictionary_t known_expr_funcs,
               const std::vector<std::string> &storeobjstr,
               const std::vector<std::string> &filesToKeep,
               const std::map<std::string, std::string> &userVariables)
    {
        this->global_rank_ = this->comm_->globalRank();

        this->parseInputFile(known_expr_funcs, false);

        MPI_Barrier(MPI_COMM_WORLD);

        // here the control flow starts to diverge
        if      ( this->comm_->isOptimizer() ) { startSampler(); }
        else if ( this->comm_->isWorker()    ) { startWorker(storeobjstr, filesToKeep, userVariables); }
        else if ( this->comm_->isPilot()     ) { this->startPilot();     }
    }

    virtual
    void startSampler() {

        std::ostringstream os;
        os << "\033[01;35m" << "  " << this->global_rank_ << " (PID: " << getpid() << ") ▶ Sampler"
           << "\e[0m" << std::endl;
        std::cout << os.str() << std::flush;

        boost::scoped_ptr<Opt_t> opt(
                                     new Opt_t(sampleMethods_m, this->objectives_, this->dvars_,
                                               this->comm_->getBundle(), this->cmd_args_));
        opt->initialize();

        std::cout << "Stop Sampler.." << std::endl;
    }

    using  Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t>::startWorker;
    void startWorker(const std::vector<std::string> &storeobjstr,
                     const std::vector<std::string> &filesToKeep,
                     const std::map<std::string, std::string> &userVariables)
    {
        std::ostringstream os;
        os << "\033[01;35m" << "  " << this->global_rank_ << " (PID: " << getpid() << ") ▶ Worker"
           << "\e[0m" << std::endl;
        std::cout << os.str() << std::flush;

        size_t pos = this->input_file_.find_last_of("/");
        std::string tmplfile = this->input_file_;
        if (pos != std::string::npos)
            tmplfile = this->input_file_.substr(pos+1);
        pos = tmplfile.find(".");
        std::string simName = tmplfile.substr(0,pos);

        boost::scoped_ptr< SampleWorker<Sim_t> > w(
                                                   new SampleWorker<Sim_t>(this->objectives_, this->constraints_, simName,
                                                                           this->comm_->getBundle(), this->cmd_args_,
                                                                           storeobjstr, filesToKeep, userVariables));

        std::cout << "Stop Worker.." << std::endl;
    }

    virtual
    void postPoll() {

        // terminating all workers is tricky since we do not know their state.
        // All workers are notified (to terminate) when opt has converged and
        // all workers are idle.
        bool all_worker_idle = true;

        // in the case where new requests became available after worker
        // delivered last results (and switched to idle state).
        for(int i = 0; i < this->total_available_workers_; i++) {

            if (i == this->my_rank_in_worker_comm_) continue;

            if (this->is_worker_idle_[i] && request_queue_.size() > 0)
                sendNewJobToWorker(i);

            all_worker_idle = all_worker_idle && this->is_worker_idle_[i];
        }

        // when all workers have been notified we can stop polling
        if (all_worker_idle && this->has_opt_converged_) {
            this->continue_polling_ = false;
            int dummy = 0;
            for(int worker = 0; worker < this->total_available_workers_; worker++) {
                MPI_Request req;
                MPI_Isend(&dummy, 1, MPI_INT, worker,
                          MPI_STOP_TAG, this->worker_comm_, &req);
            }
        }
    }


    virtual
    void sendNewJobToWorker(int worker) /*override*/ {

        // no new jobs once our opt has converged
        if (this->has_opt_converged_) return;

        JobIter_t job = request_queue_.begin();
        size_t jid = job->first;

        Param_t job_params = job->second;
        MPI_Send(&jid, 1, MPI_UNSIGNED_LONG, worker, MPI_WORK_JOBID_TAG, this->worker_comm_);
        MPI_Send_params(job_params, worker, this->worker_comm_);

        running_job_list_.insert(std::pair<size_t,
                                 Param_t >(job->first, job->second));
        request_queue_.erase(jid);
        this->is_worker_idle_[worker] = false;

        std::ostringstream dump;
        dump << "sent job with ID " << jid << " to worker " << worker
             << std::endl;
        this->job_trace_->log(dump);

    }


    virtual
    bool onMessage(MPI_Status status, size_t recv_value) /*override*/ {

        MPITag_t tag = MPITag_t(status.MPI_TAG);
        switch(tag) {

        case WORKER_FINISHED_TAG: {

            size_t job_id = recv_value;

            size_t dummy = 1;
            MPI_Send(&dummy, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE,
                     MPI_WORKER_FINISHED_ACK_TAG, this->worker_comm_);

            reqVarContainer_t res;
            MPI_Recv_reqvars(res, status.MPI_SOURCE, this->worker_comm_);

            running_job_list_.erase(job_id);
            this->is_worker_idle_[status.MPI_SOURCE] = true;

            std::ostringstream dump;
            dump << "worker finished job with ID " << job_id << std::endl;
            this->job_trace_->log(dump);


            // sampler already terminated, cannot accept new messages
            if (this->has_opt_converged_) return true;

            int opt_master_rank = this->comm_->getLeader();
            MPI_Send(&job_id, 1, MPI_UNSIGNED_LONG, opt_master_rank,
                     MPI_OPT_JOB_FINISHED_TAG, this->opt_comm_);

            MPI_Send_reqvars(res, opt_master_rank, this->opt_comm_);

            // we keep worker busy _after_ results have been sent to sampler
            if (request_queue_.size() > 0)
                sendNewJobToWorker(status.MPI_SOURCE);

            return true;
        }

        case OPT_NEW_JOB_TAG: {

            size_t job_id = recv_value;
            int opt_master_rank = this->comm_->getLeader();

            Param_t job_params;
            MPI_Recv_params(job_params, (size_t)opt_master_rank, this->opt_comm_);

            request_queue_.insert(
                                  std::pair<size_t, Param_t >(
                                                              job_id, job_params));

            std::ostringstream dump;
            dump << "new opt job with ID " << job_id << std::endl;
            this->job_trace_->log(dump);

            return true;
        }

        case OPT_CONVERGED_TAG: {
            return this->stop();
        }

        case WORKER_STATUSUPDATE_TAG: {
            this->is_worker_idle_[status.MPI_SOURCE] = true;
            return true;
        }

        default: {
            std::string msg = "(Pilot) Error: unexpected MPI_TAG: ";
            msg += status.MPI_TAG;
            throw OptPilotException("SamplePilot::onMessage", msg);
        }
        }
    }

private:
    std::map< std::string,
              std::shared_ptr<SamplingMethod>
              > sampleMethods_m;
};

#endif
