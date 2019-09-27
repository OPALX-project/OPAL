#ifndef __PILOT_H__
#define __PILOT_H__

#include <mpi.h>
#include <iostream>
#include <string>

#include "boost/smart_ptr.hpp"
//#include "boost/dynamic_bitset.hpp"

#include "Comm/MasterNode.h"
#include "Comm/CommSplitter.h"

#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Util/InputFileParser.h"
#include "Util/OptPilotException.h"

#include "Pilot/Poller.h"
#include "Pilot/Worker.h"
#include "Optimizer/Optimizer.h"

#include "Util/Trace/Trace.h"
#include "Util/Trace/FileSink.h"
#include "Util/Trace/TraceComponent.h"

#include "Expression/Parser/function.hpp"


/*! \mainpage Opt-Pilot
 *
 * \section intro_sec Introduction
 *
 * This is a general-purpose framework for simulation-based multi-objective
 * optimization methods that allows the automatic investigation of Pareto
 * fronts.
 *
 * The implementation is based on a Master/Slave paradigm, employing several
 * masters and groups of workers. To tackle the emerging huge problems
 * efficiently, we employ network topology-aware mappings of masters and
 * slaves.
 *
 * The framework is easy to extend and use with other simulation-based forward
 * solvers, optimization algorithms and network mapping strategies.
 * Currently the code contains bindings for:
 *
 *   - OPAL forward simulations
 *   - EA based optimizer (see PISA) using a NSGA-II selector
 *   - one master
 *   - multiple masters
 *   - a number of co-workers for running the simulations
 *
 * See the README file for build instructions.
*/



/**
 *  \class Pilot
 *  \brief The Optimization Pilot (Master): Coordinates requests by optimizer
 *  to workers and reports results back on a given communicator.
 *
 *  Every worker thread notifies the master here if idle or not. When
 *  available the master dispatches one of the pending simulations to the
 *  worker who will run the specified simulation and report results back to
 *  the master. The Optimizer class will poll the scheduler to check if some
 *  (or all) results are available and continue to optimize and request new
 *  simulation results.
 *
 *  @see Worker
 *  @see Optimizer
 *
 *  @tparam Input_t type of the input file parser
 *  @tparam Opt_t type of the optimizer
 *  @tparam Sim_t type of the simulation
 *  @tparam SolPropagationGraph_t strategy to distribute solution between
 *          master islands
 *  @tparam Comm_t comm splitter strategy
 */
template <
          class Input_t
        , class Opt_t
        , class Sim_t
        , class SolPropagationGraph_t
        , class Comm_t
>
class Pilot : protected Poller {

public:

    // constructor only for Pilot classes inherited from this class
    // they have their own setup function
    Pilot(CmdArguments_t args, boost::shared_ptr<Comm_t> comm,
          const DVarContainer_t &dvar)
        : Poller(comm->mpiComm())
        , comm_(comm)
        , cmd_args_(args)
        , dvars_(dvar)
    {
        // do nothing
    }

    Pilot(CmdArguments_t args, boost::shared_ptr<Comm_t> comm,
          functionDictionary_t known_expr_funcs,
          std::vector<std::string> filesToKeep = {})
        : Poller(comm->mpiComm())
        , comm_(comm)
        , cmd_args_(args)
    {
        setup(known_expr_funcs, filesToKeep);
    }

    Pilot(CmdArguments_t args, boost::shared_ptr<Comm_t> comm,
          functionDictionary_t known_expr_funcs,
          const DVarContainer_t &dvar,
          const Expressions::Named_t &obj,
          const Expressions::Named_t &cons,
          std::vector<std::string> filesToKeep = {})
        : Poller(comm->mpiComm())
        , comm_(comm)
        , cmd_args_(args)
        , objectives_(obj)
        , constraints_(cons)
        , dvars_(dvar)
    {
        setup(known_expr_funcs, filesToKeep);
    }

    ~Pilot()
    {}


protected:

    /// MPI communicator used for messages to/from worker
    MPI_Comm worker_comm_;
    /// MPI communicator used for messages to/from optimizer
    MPI_Comm opt_comm_;
    /// MPI communicator used for messages between all pilots
    MPI_Comm coworker_comm_;

    boost::shared_ptr<Comm_t> comm_;
    CmdArguments_t cmd_args_;

    int global_rank_;
    int my_rank_in_worker_comm_;
    int my_rank_in_opt_comm_;

    int num_coworkers_;

    typedef MasterNode< typename Opt_t::SolutionState_t,
                        SolPropagationGraph_t > MasterNode_t;
    boost::scoped_ptr< MasterNode_t > master_node_;

    /// input file for simulation with embedded optimization problem
    std::string input_file_;

    int  total_available_workers_;
    bool has_opt_converged_;
    bool continue_polling_;

    Expressions::Named_t objectives_;
    Expressions::Named_t constraints_;
    DVarContainer_t      dvars_;

    // keep track of state of all workers
    std::vector<bool> is_worker_idle_;
    //boost::dynamic_bitset<> is_worker_idle_;

    /// keep track of requests and running jobs
    typedef std::map<size_t, std::pair<Param_t, reqVarContainer_t> > Jobs_t;
    typedef Jobs_t::iterator JobIter_t;
    Jobs_t  running_job_list_;
    Jobs_t  request_queue_;

    //DEBUG
    boost::scoped_ptr<Trace> job_trace_;

private:
    void setup(functionDictionary_t known_expr_funcs,
               const std::vector<std::string> &filesToKeep)
    {
        global_rank_ = comm_->globalRank();

        if(global_rank_ == 0) {
            std::cout << "\033[01;35m";
            std::cout << "             _                _ _       _   " << std::endl;
            std::cout << "            | |              (_) |     | |  " << std::endl;
            std::cout << "  ___  _ __ | |_ ______ _ __  _| | ___ | |_ " << std::endl;
            std::cout << " / _ \\| '_ \\| __|______| '_ \\| | |/ _ \\| __|" << std::endl;
            std::cout << "| (_) | |_) | |_       | |_) | | | (_) | |_ " << std::endl;
            std::cout << " \\___/| .__/ \\__|      | .__/|_|_|\\___/ \\__|" << std::endl;
            std::cout << "      | |              | |                  " << std::endl;
            std::cout << "      |_|              |_|                  " << std::endl;
	    // ADA            std::cout << "☷ Version: \t"    << PACKAGE_VERSION << std::endl;
            //std::cout << "☷ Git: \t\t"      << GIT_VERSION     << std::endl;
            //std::cout << "☷ Build Date: \t" << BUILD_DATE      << std::endl;
            std::cout << "\e[0m";
            std::cout << std::endl;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        parseInputFile(known_expr_funcs);

        // here the control flow starts to diverge
        if      ( comm_->isOptimizer() ) { startOptimizer(); }
        else if ( comm_->isWorker()    ) { startWorker(filesToKeep);    }
        else if ( comm_->isPilot()     ) { startPilot();     }
    }

protected:

    void parseInputFile(functionDictionary_t known_expr_funcs) {

        try {
            input_file_ = cmd_args_->getArg<std::string>("inputfile", true);
        } catch (OptPilotException &e) {
            std::cout << "Could not find 'inputfile' in arguments.. Aborting."
                << std::endl;
            MPI_Abort(comm_m, -101);
        }

        if(objectives_.size() == 0 || dvars_.size() == 0) {
            throw OptPilotException("Pilot::Pilot()",
                    "No objectives or dvars specified");
        }

        if(global_rank_ == 0) {
            std::ostringstream os;
            os << "\033[01;35m";
            os << "  ✔ " << objectives_.size()
               << " objectives" << std::endl;
            os << "  ✔ " << constraints_.size()
               << " constraints" << std::endl;
            os << "  ✔ " << dvars_.size()
               << " dvars" << std::endl;
            os << "\e[0m";
            os << std::endl;
            std::cout << os.str() << std::flush;
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    virtual
    void startOptimizer() {

        std::ostringstream os;
        os << "\033[01;35m" << "  " << global_rank_ << " ▶ Opt"
           << "\e[0m" << std::endl;
        std::cout << os.str() << std::flush;

        boost::scoped_ptr<Opt_t> opt(
                new Opt_t(objectives_, constraints_, dvars_, objectives_.size(),
                    comm_->getBundle(), cmd_args_));
        opt->initialize();

        std::cout << "Stop Opt.." << std::endl;
    }

    virtual
    void startWorker(const std::vector<std::string> &filesToKeep) {

        std::ostringstream os;
        os << "\033[01;35m" << "  " << global_rank_ << " ▶ Worker"
           << "\e[0m" << std::endl;
        std::cout << os.str() << std::flush;

        size_t pos = input_file_.find_last_of("/");
        std::string tmplfile = input_file_;
        if(pos != std::string::npos)
            tmplfile = input_file_.substr(pos+1);
        pos = tmplfile.find(".");
        std::string simName = tmplfile.substr(0,pos);

        boost::scoped_ptr< Worker<Sim_t> > w(
                new Worker<Sim_t>(objectives_, constraints_, simName,
                    comm_->getBundle(), cmd_args_, filesToKeep));

        std::cout << "Stop Worker.." << std::endl;
    }

    virtual
    void startPilot() {

        std::ostringstream os;
        os << "\033[01;35m" << "  " << global_rank_ << " ▶ Pilot"
           << "\e[0m" << std::endl;
        std::cout << os.str() << std::flush;

        // Traces
        std::ostringstream trace_filename;
        trace_filename << "pilot.trace." << comm_->getBundle().island_id;
        job_trace_.reset(new Trace("Optimizer Job Trace"));
        job_trace_->registerComponent( "sink",
            boost::shared_ptr<TraceComponent>(new FileSink(trace_filename.str())));

        worker_comm_ = comm_->getBundle().worker;
        opt_comm_ = comm_->getBundle().opt;
        coworker_comm_ = comm_->getBundle().world;

        my_rank_in_worker_comm_ = 0;
        MPI_Comm_rank(worker_comm_, &my_rank_in_worker_comm_);
        my_rank_in_opt_comm_ = 0;
        MPI_Comm_rank(opt_comm_, &my_rank_in_opt_comm_);

        total_available_workers_ = 0;
        MPI_Comm_size(worker_comm_, &total_available_workers_);
        is_worker_idle_.resize(total_available_workers_, true);
        is_worker_idle_[my_rank_in_worker_comm_] = false;

        // setup master network
        num_coworkers_ = 0;
        MPI_Comm_size(coworker_comm_, &num_coworkers_);
        if(num_coworkers_ > 1) {
            //FIXME: proper upper bound for window size
            int alpha = cmd_args_->getArg<int>("initialPopulation", false);
            int opt_size = objectives_.size() + constraints_.size();
            int overhead = 10;
            size_t upperbound_buffer_size =
                sizeof(double) * alpha * (1 + opt_size) * 1000
                + overhead;
            master_node_.reset(
                new MasterNode< typename Opt_t::SolutionState_t,
                                SolPropagationGraph_t >(
                    coworker_comm_, upperbound_buffer_size, objectives_.size(),
                    comm_->getBundle().island_id));
        }

        has_opt_converged_ = false;
        continue_polling_  = true;
        run();

        std::cout << "Stop Pilot.." << std::endl;
    }

    virtual
    void setupPoll()
    {}

    virtual
    void prePoll()
    {}

    virtual
    void onStop()
    {}

    virtual
    void postPoll() {

        // terminating all workers is tricky since we do not know their state.
        // All workers are notified (to terminate) when opt has converged and
        // all workers are idle.
        bool all_worker_idle = true;

        // in the case where new requests became available after worker
        // delivered last results (and switched to idle state).
        for(int i = 0; i < total_available_workers_; i++) {

            if(i == my_rank_in_worker_comm_) continue;

            all_worker_idle = all_worker_idle && is_worker_idle_[i];

            if(is_worker_idle_[i] && request_queue_.size() > 0)
                sendNewJobToWorker(i);
        }

        // when all workers have been notified we can stop polling
        if(all_worker_idle && has_opt_converged_) {
            continue_polling_ = false;
            int dummy = 0;
            for(int worker = 0; worker < total_available_workers_; worker++) {
                MPI_Request req;
                MPI_Isend(&dummy, 1, MPI_INT, worker,
                          MPI_STOP_TAG, worker_comm_, &req);
            }
        }
    }


    virtual
    void sendNewJobToWorker(int worker) {

        // no new jobs once our opt has converged
        if(has_opt_converged_) return;

        JobIter_t job = request_queue_.begin();
        size_t jid = job->first;

        Param_t job_params = job->second.first;
        MPI_Send(&jid, 1, MPI_UNSIGNED_LONG, worker, MPI_WORK_JOBID_TAG, worker_comm_);
        MPI_Send_params(job_params, worker, worker_comm_);

        //reqVarContainer_t job_reqvars = job->second.second;
        //MPI_Send_reqvars(job_reqvars, worker, worker_comm_);

        running_job_list_.insert(std::pair<size_t,
                std::pair<Param_t, reqVarContainer_t> >(job->first, job->second));
        request_queue_.erase(jid);
        is_worker_idle_[worker] = false;

        std::ostringstream dump;
        dump << "sent job with ID " << jid << " to worker " << worker
             << std::endl;
        job_trace_->log(dump);

    }


    virtual
    bool onMessage(MPI_Status status, size_t recv_value){

        MPITag_t tag = MPITag_t(status.MPI_TAG);
        switch(tag) {

        case WORKER_FINISHED_TAG: {

            size_t job_id = recv_value;

            size_t dummy = 1;
            MPI_Send(&dummy, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE,
                     MPI_WORKER_FINISHED_ACK_TAG, worker_comm_);

            MPI_Status st;
            MPI_Probe(status.MPI_SOURCE, MPI_WORKER_DIRECTORY_TAG, worker_comm_, &st);
            int len = 0;
            MPI_Get_count(&st, MPI_CHAR, &len);
            std::unique_ptr<char> buf(new char[len]);
            MPI_Recv(buf.get(), len, MPI_CHAR, status.MPI_SOURCE,
                     MPI_WORKER_DIRECTORY_TAG, worker_comm_, &st);
            std::string simdir(buf.get(), len);


            reqVarContainer_t res;
            MPI_Recv_reqvars(res, status.MPI_SOURCE, worker_comm_);

            running_job_list_.erase(job_id);
            is_worker_idle_[status.MPI_SOURCE] = true;

            std::ostringstream dump;
            dump << "worker finished job with ID " << job_id << std::endl;
            job_trace_->log(dump);


            // optimizer already terminated, cannot accept new messages
            if(has_opt_converged_) return true;

            int opt_master_rank = comm_->getLeader();
            MPI_Send(&job_id, 1, MPI_UNSIGNED_LONG, opt_master_rank,
                     MPI_OPT_JOB_FINISHED_TAG, opt_comm_);

            MPI_Send(simdir.c_str(), simdir.size(), MPI_CHAR, opt_master_rank,
                     MPI_OPT_DIRECTORY_TAG, opt_comm_);

            MPI_Send_reqvars(res, opt_master_rank, opt_comm_);

            // we keep worker busy _after_ results have been sent to optimizer
            if(request_queue_.size() > 0)
                sendNewJobToWorker(status.MPI_SOURCE);

            return true;
        }

        case OPT_NEW_JOB_TAG: {

            size_t job_id = recv_value;
            int opt_master_rank = comm_->getLeader();

            Param_t job_params;
            MPI_Recv_params(job_params, (size_t)opt_master_rank, opt_comm_);

            reqVarContainer_t reqVars;
            //MPI_Recv_reqvars(reqVars, (size_t)opt_master_rank, job_size, opt_comm_);

            std::pair<Param_t, reqVarContainer_t> job =
                std::pair<Param_t, reqVarContainer_t>(job_params, reqVars);
            request_queue_.insert(
                                  std::pair<size_t, std::pair<Param_t, reqVarContainer_t> >(
                                                                                            job_id, job));

            std::ostringstream dump;
            dump << "new opt job with ID " << job_id << std::endl;
            job_trace_->log(dump);

            return true;
        }

        case EXCHANGE_SOL_STATE_TAG: {

            if(num_coworkers_ <= 1) return true;

            std::ostringstream dump;
            dump << "starting solution exchange.. " << status.MPI_SOURCE << std::endl;
            job_trace_->log(dump);

            // we start by storing or local solution state
            size_t buffer_size = recv_value;
            int opt_master_rank = status.MPI_SOURCE; //comm_->getLeader();

            char *buffer = new char[buffer_size];
            MPI_Recv(buffer, buffer_size, MPI_CHAR, opt_master_rank,
                     MPI_EXCHANGE_SOL_STATE_DATA_TAG, opt_comm_, &status);
            master_node_->store(buffer, buffer_size);
            delete[] buffer;

            dump.clear();
            dump.str(std::string());
            dump << "getting " << buffer_size << " bytes from OPT "
                 << opt_master_rank << std::endl;
            job_trace_->log(dump);

            // and then continue collecting all other solution states
            std::ostringstream states;
            master_node_->collect(states);
            buffer_size = states.str().length();

            dump.clear();
            dump.str(std::string());
            dump << "collected solution states of other PILOTS: "
                 << buffer_size << " bytes" << std::endl;
            job_trace_->log(dump);

            // send collected solution states to optimizer;
            MPI_Send(&buffer_size, 1, MPI_UNSIGNED_LONG, opt_master_rank,
                     MPI_EXCHANGE_SOL_STATE_RES_SIZE_TAG, opt_comm_);

            buffer = new char[buffer_size];
            memcpy(buffer, states.str().c_str(), buffer_size);
            MPI_Send(buffer, buffer_size, MPI_CHAR, opt_master_rank,
                     MPI_EXCHANGE_SOL_STATE_RES_TAG, opt_comm_);

            dump.clear();
            dump.str(std::string());
            dump << "sent set of new solutions to OPT" << std::endl;
            job_trace_->log(dump);

            delete[] buffer;

            return true;
        }

        case OPT_CONVERGED_TAG: {
            return stop();
        }

        case WORKER_STATUSUPDATE_TAG: {
            is_worker_idle_[status.MPI_SOURCE] = true;
            return true;
        }

        default: {
            std::string msg = "(Pilot) Error: unexpected MPI_TAG: ";
            msg += status.MPI_TAG;
            throw OptPilotException("Pilot::onMessage", msg);
        }
        }
    }

    bool stop(bool isOpt = true) {

        if(has_opt_converged_) return true;

        has_opt_converged_ = true;
        request_queue_.clear();
        size_t dummy = 0;
        MPI_Request req;
        MPI_Isend(&dummy, 1, MPI_UNSIGNED_LONG, comm_->getLeader(), MPI_STOP_TAG, opt_comm_, &req);

        if(! isOpt) return true;
        if(num_coworkers_ <= 1) return true;

        if(! cmd_args_->getArg<bool>("one-pilot-converge", false, false))
            return true;

        // propagate converged message to other pilots
        // FIXME what happens if two island converge at the same time?
        int my_rank = 0;
        MPI_Comm_rank(coworker_comm_, &my_rank);
        for(int i=0; i < num_coworkers_; i++) {
            if(i == my_rank) continue;
            MPI_Request req;
            MPI_Isend(&dummy, 1, MPI_UNSIGNED_LONG, i, OPT_CONVERGED_TAG, coworker_comm_, &req);
        }

        return true;
    }


    // we overwrite run here to handle polling on two different communicators
    //XXX: would be nice to give the poller interface an array of comms and
    //     listeners to be called..
    void run() {

        MPI_Request opt_request;
        MPI_Request worker_request;
        MPI_Status status;
        int flag = 0;
        size_t recv_value_worker = 0;
        size_t recv_value_opt = 0;

        setupPoll();

        MPI_Irecv(&recv_value_opt, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE,
                  MPI_ANY_TAG, opt_comm_, &opt_request);
        MPI_Irecv(&recv_value_worker, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE,
                  MPI_ANY_TAG, worker_comm_, &worker_request);

        bool pending_opt_request    = true;
        bool pending_worker_request = true;
        bool pending_pilot_request  = false;

        MPI_Request pilot_request;
        size_t recv_value_pilot = 0;
        if(cmd_args_->getArg<bool>("one-pilot-converge", false, false)) {
            MPI_Irecv(&recv_value_pilot, 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE,
                      MPI_ANY_TAG, coworker_comm_, &pilot_request);
            pending_pilot_request = true;
        }

        while(continue_polling_) {

            prePoll();

            if(opt_request != MPI_REQUEST_NULL) {
                MPI_Test(&opt_request, &flag, &status);
                if(flag) {
                    pending_opt_request = false;
                    if(status.MPI_TAG == MPI_STOP_TAG) {
                        return;
                    } else {
                        if(onMessage(status, recv_value_opt)) {
                            MPI_Irecv(&recv_value_opt, 1, MPI_UNSIGNED_LONG,
                                      MPI_ANY_SOURCE, MPI_ANY_TAG, opt_comm_,
                                      &opt_request);
                            pending_opt_request = true;
                        } else
                            return;
                    }
                }
            }

            if(worker_request != MPI_REQUEST_NULL) {
                MPI_Test(&worker_request, &flag, &status);
                if(flag) {
                    pending_worker_request = false;
                    if(status.MPI_TAG == MPI_STOP_TAG) {
                        return;
                    } else {
                        if(onMessage(status, recv_value_worker)) {
                            MPI_Irecv(&recv_value_worker, 1,
                                      MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                                      worker_comm_, &worker_request);
                            pending_worker_request = true;
                        } else
                            return;
                    }
                }
            }

            if(cmd_args_->getArg<bool>("one-pilot-converge", false, false)) {
                if(pilot_request != MPI_REQUEST_NULL) {
                    MPI_Test(&pilot_request, &flag, &status);
                    if(flag) {
                        pending_pilot_request = false;
                        if(status.MPI_TAG == OPT_CONVERGED_TAG) {
                            stop(false);
                        } else {
                            MPI_Irecv(&recv_value_pilot, 1,
                                      MPI_UNSIGNED_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG,
                                      coworker_comm_, &pilot_request);
                            pending_pilot_request = true;
                        }
                    }
                }
            }

            postPoll();
        }

        if(pending_opt_request)     MPI_Cancel( &opt_request );
        if(pending_worker_request)  MPI_Cancel( &worker_request );
        if(pending_pilot_request)   MPI_Cancel( &pilot_request );
    }

};

#endif