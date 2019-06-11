#include <iostream>
#include <sstream>
#include <cmath>
#include <cstring>
#include <set>
#include <cstdlib>
#include <string>
#include <limits>

#include <sys/stat.h>

#include "boost/algorithm/string.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>

#include "OPALconfig.h"
#include "Utilities/Util.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "Util/OptPilotException.h"
#include "Util/MPIHelper.h"

#include "Util/Trace/TraceComponent.h"
#include "Util/Trace/Timestamp.h"
#include "Util/Trace/FileSink.h"


template< template <class> class CO, template <class> class MO >
FixedPisaNsga2<CO, MO>::FixedPisaNsga2(
                           Expressions::Named_t objectives,
                           Expressions::Named_t constraints,
                           DVarContainer_t dvars,
                           size_t dim, Comm::Bundle_t comms,
                           CmdArguments_t args,
                           std::vector<double> hypervolRef,
                           int nrWorkerGroups)
             : Optimizer(comms.opt)
             , statistics_(new Statistics<size_t>("individuals"))
             , comms_(comms)
             , objectives_m(objectives)
             , constraints_m(constraints)
             , dvars_m(dvars)
             , args_m(args)
             , dim_m(dim)
             , num_workergroups_m(nrWorkerGroups)
             , hvol_ref_m(hypervolRef)
{
    //FIXME: proper rand gen initialization (use boost?!)
    srand(time(NULL) + comms_.island_id);

    dump_freq_m      = args->getArg<int>("dump-freq", 1, false);
    dump_offspring_m = args->getArg<bool>("dump-offspring", true, false);
    dump_dat_m       = args->getArg<bool>("dump-dat", false, false);

    maxGenerations_m = args->getArg<int>("maxGenerations", true);
    resultFile_m     = args->getArg<std::string>("outfile", "-th_generation.dat", false);
    resultDir_m      = args->getArg<std::string>("outdir", "generations", false);

    // create output directory if it does not exists
    struct stat dirInfo;
    if(stat(resultDir_m.c_str(),&dirInfo) != 0)
        mkdir((const char*)(resultDir_m.c_str()), 0777);

    // solution exchange frequency
    exchangeSolStateFreq_m = args->getArg<size_t>("sol-synch", 0, false);

    // convergence arguments
    hvol_eps_           = args->getArg<double>("epsilon", 1e-3, false);
    expected_hvol_      = args->getArg<double>("expected-hypervol", 0.0,
                                               false);
    conv_hvol_progress_ = args->getArg<double>("conv-hvol-prog", 0.0, false);
    current_hvol_       = std::numeric_limits<double>::max();
    hvol_progress_      = std::numeric_limits<double>::max();

    lambda_m              = args->getArg<int>("num-ind-gen", 2, false);
    alpha_m               = args->getArg<int>("initialPopulation", true);
    //mu_m                = lambda_m;

    // initial population arguments
    file_start_m          = args_m->getArg<std::string>("start-population",     "",    false);
    initialOptimization_m = args_m->getArg<bool>       ("initial-optimization", false, false);
    birthControl_m        = args_m->getArg<bool>       ("birth-control",        false, false);

    file_param_descr_ = "%ID,";

    Expressions::Named_t::iterator it;
    for(it = objectives_m.begin(); it != objectives_m.end(); it++)
        file_param_descr_ += '%' + it->first + ',';

    file_param_descr_ += " DVAR: ";

    DVarContainer_t::iterator itr;
    std::vector<std::string> dNames;

    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        std::string dName = boost::get<VAR_NAME>(itr->second);
        file_param_descr_ += '%' + dName + ',';
        dNames.push_back(dName);
        dVarBounds_m.push_back(
                               std::pair<double, double>
                               (boost::get<LOWER_BOUND>(itr->second),
                                boost::get<UPPER_BOUND>(itr->second)));
    }
    file_param_descr_ = file_param_descr_.substr(0,
                                                 file_param_descr_.size()-1);

    // setup variator
    variator_m.   reset(new Variator_t(dNames, dVarBounds_m, constraints_m, args));
    paretoFront_m.reset(new Population_t());

    // Traces and statistics
    std::ostringstream trace_filename;
    trace_filename << "opt.trace." << comms_.island_id;
    job_trace_.reset(new Trace("Optimizer Job Trace"));
    job_trace_->registerComponent( "sink",
            boost::shared_ptr<TraceComponent>(
                new FileSink(trace_filename.str())));

    std::ostringstream prog_filename;
    prog_filename << "opt.progress." << comms_.island_id;
    progress_.reset(new Trace("Optimizer Progress"));
    progress_->registerComponent( "timestamp",
            boost::shared_ptr<TraceComponent>(new Timestamp()));
    progress_->registerComponent( "sink",
            boost::shared_ptr<TraceComponent>(
                new FileSink(prog_filename.str())));

    statistics_->registerStatistic("accepted", 0);
    statistics_->registerStatistic("infeasible", 0);
}

template<
      template <class> class CO
    , template <class> class MO
>
FixedPisaNsga2<CO, MO>::~FixedPisaNsga2()
{}

template<
      template <class> class CO
    , template <class> class MO
>
void FixedPisaNsga2<CO, MO>::initialize() {

    curState_m = Initialize;
    initialized_m = false;

    // start poll loop
    run_clock_start_ = boost::chrono::system_clock::now();
    last_clock_      = boost::chrono::system_clock::now();
    run();

    // run has ended
    bool compHyvol = (objectives_m.size() > (hyper_opt / 2 + 1));
    if (compHyvol)
        current_hvol_ =
          variator_m->population()->computeHypervolume(comms_.island_id, hvol_ref_m);

    boost::chrono::duration<double> total =
        boost::chrono::system_clock::now() - run_clock_start_;
    std::ostringstream stats;
    stats << "__________________________________________" << std::endl;
    stats << "GENERATION " <<  act_gen << std::endl;
    stats << "TOTAL = " << total.count() << "s" << std::endl;
    if (compHyvol)
        stats << "HYPERVOLUME = " << current_hvol_ << std::endl;
    stats << "time per accepted ind   = "
          << total.count()/statistics_->getStatisticValue("accepted")
          << std::endl;
    stats << "time per infeasible ind = "
          << (statistics_->getStatisticValue("infeasible") ? total.count()/statistics_->getStatisticValue("infeasible") : 0)
          << std::endl;
    statistics_->dumpStatistics(stats);
    stats << "__________________________________________" << std::endl;
    progress_->log(stats);
}

template< template <class> class CO, template <class> class MO >
bool FixedPisaNsga2<CO, MO>::onMessage(MPI_Status status, size_t length) {

    MPITag_t tag = MPITag_t(status.MPI_TAG);
    switch(tag) {

    case EXCHANGE_SOL_STATE_RES_SIZE_TAG: {

        size_t buf_size = length;
        size_t pilot_rank = status.MPI_SOURCE;

        std::ostringstream dump;
        dump << "new results from other cores " << buf_size << std::endl;
        job_trace_->log(dump);

        char *buffer = new char[buf_size];
        MPI_Recv(buffer, buf_size, MPI_CHAR, pilot_rank,
                 MPI_EXCHANGE_SOL_STATE_RES_TAG, comms_.opt, &status);

        dump.clear();
        dump.str(std::string());
        dump << "got results from other cores " << buf_size << std::endl;
        job_trace_->log(dump);

        SolutionState_t new_states;
        std::istringstream is(buffer);
        boost::archive::text_iarchive ia(is);
        ia >> new_states;
        delete[] buffer;

        std::set<unsigned int> new_state_ids;
        for (Individual_t ind : new_states) {

            // only insert individual if not already in population
            if(variator_m->population()->isRepresentedInPopulation(ind.genes_m))
                continue;

            individual new_ind(new Individual_t);
            new_ind->genes_m = ind.genes_m;
            new_ind->objectives_m = ind.objectives_m;

            //XXX:   can we pass more than lambda_m files to selector?
            unsigned int id =
                variator_m->population()->add_individual(new_ind);
            finishedBuffer_m.push_back(id);

            dump.clear();
            dump.str(std::string());
            dump << "individual (ID: " << id
                 << ") successfully migrated from another island" << std::endl;
            job_trace_->log(dump);
        }

        return true;
    }

    case REQUEST_FINISHED: {

        unsigned int jid = static_cast<unsigned int>(length);
        typename std::map<size_t, individual >::iterator it;
        it = jobmapping_m.find(jid);

        std::ostringstream dump;
        dump << "job with ID " << jid << " delivered results" << std::endl;
        job_trace_->log(dump);

        if(it == jobmapping_m.end()) {
            dump << "\t |-> NOT FOUND!" << std::endl;
            job_trace_->log(dump);
            std::cout << "NON-EXISTING JOB with ID = " << jid << std::endl;
            throw OptPilotException("FixedPisaNsga2::onMessage",
                                    "non-existing job");
        }

        individual ind = it->second;
        jobmapping_m.erase(it);

        //size_t dummy = 1;
        //MPI_Send(&dummy, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE,
        //         MPI_WORKER_FINISHED_ACK_TAG, comms_.listen);

        reqVarContainer_t res;
        MPI_Recv_reqvars(res, status.MPI_SOURCE, comms_.opt);

        ind->objectives_m.clear();

        //XXX: check order of genes
        reqVarContainer_t::iterator itr = res.begin();
        for(; itr != res.end(); ++itr) {
            // mark invalid if expression could not be evaluated or constraint does not hold
            if(!itr->second.is_valid || (itr->second.value.size() > 1 && !itr->second.value[0])) {
                std::ostringstream dump;
                if (!itr->second.is_valid) {
                    dump << "invalid individual, objective or constraint \"" << itr->first
                         << "\" failed to be evaluated correctly"
                         << std::endl;
                } else {
                    dump << "invalid individual, constraint \"" << itr->first
                         << "\" failed to yield true; result: " << itr->second.value[1]
                         << std::endl;
                }
                job_trace_->log(dump);
                variator_m->infeasible(ind);
                statistics_->changeStatisticBy("infeasible", 1);
                dispatch_forward_solves();
                return true;
            } else {
                 // update objective value for valid objective
                 if(itr->second.value.size() == 1)
                    ind->objectives_m.push_back(itr->second.value[0]);
            }
        }

        finishedBuffer_m.push_back(jid);
        statistics_->changeStatisticBy("accepted", 1);
        // check for pareto front
        checkParetoFront(jid);

        return true;
    }

    default: {
        std::cout << "(FixedPisaNsga2) Error: unexpected MPI_TAG: "
                  << status.MPI_TAG << std::endl;
        return false;
    }
    }
}


template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::postPoll() {

    // whenever lambda children are ready and we are in the variator phase
    // run the selector
    // std::ostringstream debug;
    // debug << "IN POST POLL: ";
    // debug << finishedBuffer_m.size() << " / " << lambda_m << std::endl;
    // debug << getStateString(curState_m) << std::endl;
    // progress_->log(debug);
    if(finishedBuffer_m.size() >= lambda_m && curState_m == Variate) {
        //std::cout << "▉";
        std::cout << "░" << std::flush;

        bool compHyvol = (objectives_m.size() > (hyper_opt / 2 + 1));
        if (compHyvol) {
            double hvol =
              variator_m->population()->computeHypervolume(comms_.island_id, hvol_ref_m);
            hvol_progress_ = fabs(current_hvol_ - hvol) / current_hvol_;
            current_hvol_ = hvol;
        }

        boost::chrono::duration<double> total =
            boost::chrono::system_clock::now() - run_clock_start_;
        boost::chrono::duration<double> dt    =
            boost::chrono::system_clock::now() - last_clock_;
        last_clock_ = boost::chrono::system_clock::now();
        std::ostringstream stats;
        stats << "__________________________________________" << std::endl;
        stats << "Arriving at generation " << act_gen + 1     << std::endl;
        stats << "dt = " << dt.count() << "s, total = " << total.count()
              << "s" << std::endl;
        if (compHyvol)
            stats << "Hypervolume = " << current_hvol_        << std::endl;
        stats << "__________________________________________" << std::endl;
        progress_->log(stats);

        // dump current generation (or their parents)
        if((act_gen + 1) % dump_freq_m == 0) {
            dumpPopulation(variator_m->population());
        }

        // we can only send all finished individuals (selector does not support
        // variable size).
        toSelectorAndCommit();

        exchangeSolutionStates();

        //XXX: at the end of the Variate state (previous runStateMachine()
        //     call) we can safely change the state.
        curState_m = Select;

        act_gen++;
    }

    // state will be reset to whatever is in the state file
    runStateMachine();
}


template< template <class> class CO , template <class> class MO >
void FixedPisaNsga2<CO, MO>::exchangeSolutionStates() {

    size_t num_masters = args_m->getArg<size_t>("num-masters", 1, false);

    if(num_masters <= 1 ||
       exchangeSolStateFreq_m == 0 ||
       act_gen % exchangeSolStateFreq_m != 0)
        return;

    int pilot_rank = comms_.master_local_pid;

    std::ostringstream os;
    boost::archive::text_oarchive oa(os);

    SolutionState_t population;
    typename std::map<unsigned int, individual >::iterator itr;
    for(itr  = variator_m->population()->begin();
        itr != variator_m->population()->end(); itr++) {

        Individual_t ind;
        ind.genes_m = std::vector<double>(itr->second->genes_m);
        ind.objectives_m = std::vector<double>(itr->second->objectives_m);
        population.push_back(ind);
    }

    oa << population;

    size_t buf_size = os.str().length();

    std::ostringstream dump;
    dump << "sending my buffer size " << buf_size << " bytes to PILOT"
         << std::endl;
    job_trace_->log(dump);

    MPI_Send(&buf_size, 1, MPI_UNSIGNED_LONG, pilot_rank,
             EXCHANGE_SOL_STATE_TAG, comms_.opt);

    char *buffer = new char[buf_size];
    memcpy(buffer, os.str().c_str(), buf_size);
    MPI_Send(buffer, buf_size, MPI_CHAR, pilot_rank,
             MPI_EXCHANGE_SOL_STATE_DATA_TAG, comms_.opt);
    delete[] buffer;

    dump.clear();
    dump.str(std::string());
    dump << "Sent " << buf_size << " bytes to PILOT" << std::endl;
    job_trace_->log(dump);
}


template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::toSelectorAndCommit() {

    to_selector_.clear();
    const size_t size = finishedBuffer_m.size();
    for(size_t i = 0; i < size; i++) {
        unsigned int id = finishedBuffer_m.front();
        to_selector_.insert(id);
        finishedBuffer_m.pop_front();
    }

    variator_m->population()->commit_individuals(to_selector_);
}


template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::dispatch_forward_solves() {

    while(variator_m->hasMoreIndividualsToEvaluate()) {

        //reqVarContainer_t reqs;
        //Expressions::Named_t::iterator it;
        //for(it = objectives_m.begin(); it != objectives_m.end(); it++) {
            //std::set<std::string> vars = it->second.getReqVars();
            //std::set<std::string>::iterator setitr;
            //for(setitr = vars.begin(); setitr != vars.end(); setitr++) {
                //if(reqs.count(*setitr) == 0) {
                    //reqVarInfo_t tmp;
                    //tmp.type = EVALUATE;
                    //tmp.value = 0.0;
                    //reqs.insert(std::pair<std::string, reqVarInfo_t>
                        //(*setitr, tmp));
                //}
            //}
        //}

        individual ind = variator_m->popIndividualToEvaluate();
        Param_t params;
        DVarContainer_t::iterator itr;
        size_t i = 0;
        for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++, i++) {
            params.insert(
                std::pair<std::string, double>
                    (boost::get<VAR_NAME>(itr->second),
                     ind->genes_m[i]));
        }

        size_t jid = static_cast<size_t>(ind->id_m);
        int pilot_rank = comms_.master_local_pid;

        // now send the request to the pilot
        MPI_Send(&jid, 1, MPI_UNSIGNED_LONG, pilot_rank, OPT_NEW_JOB_TAG, comms_.opt);

        MPI_Send_params(params, pilot_rank, comms_.opt);

        //MPI_Send_reqvars(reqs, pilot_rank, comms_.opt);

        jobmapping_m.insert(
                std::pair<size_t, individual >(jid, ind));

        std::ostringstream dump;
        dump << "dispatched simulation with ID " << jid << std::endl;
        job_trace_->log(dump);
    }
}


template<
      template <class> class CO
    , template <class> class MO
>
void FixedPisaNsga2<CO, MO>::runStateMachine() {

    switch(curState_m) {

    // State 2 of the FSM: we read mu parents
    // and create lambda children by variation of the individuals specified
    // in sel file.
    case Variate: {

        if( (maxGenerations_m > 0 && act_gen > maxGenerations_m) ||
            (hvol_progress_ < conv_hvol_progress_) ||
            (expected_hvol_ > 0.0 && fabs(current_hvol_ - expected_hvol_) < hvol_eps_)
          ) {
            curState_m = Stop;
        } else {
            if(parent_queue_.size() > 0) {
                std::vector<unsigned int> parents(parent_queue_.begin(),
                                                  parent_queue_.end());
                parent_queue_.clear();

                // remove non-surviving individuals
                //std::set<unsigned int> survivors(archive_.begin(),
                                                 //archive_.end());
                std::set<unsigned int> survivors(pp_all.begin(),
                                                 pp_all.end());

                variator_m->population()->keepSurvivors(survivors);

                // only enqueue new individuals to pilot, the poll loop will
                // feed results back to the selector.
                //FIXME: variate works on staging set!!!
                variator_m->variate(parents);
                dispatch_forward_solves();

                curState_m = Variate;
            }
        }
        break;
    }

    // State 0 of the FSM: generate initial population and write
    // information about initial population to ini file.
    case Initialize: {
        if(initialized_m == false) {
            variator_m->initial_population(alpha_m, file_start_m);
            if (initialOptimization_m == true && birthControl_m == false && file_start_m.empty()) {
                // increase population such that workers keep busy until full initial population is found
                variator_m->initial_population(num_workergroups_m, "");
            }
            dispatch_forward_solves();
            initialized_m = true;
        }

        //XXX: wait till the first alpha_m individuals are ready then start
        //     the selector
        if(finishedBuffer_m.size() >= alpha_m) {
            if(dump_offspring_m == true) { // dump first generation - no parents yet
                dumpPopulation(variator_m->population());
            }
            act_gen = 2;
            toSelectorAndCommit();
            curState_m = InitializeSelector;

            // add additional worker jobs to the queue so that single job won't keep system idle
            if (initialOptimization_m == false && birthControl_m == false && file_start_m.empty()) {
                // increase population to have all workers busy
                variator_m->initial_population(num_workergroups_m, "");
                dispatch_forward_solves();
            }
        }
        break;
    }

    // State 4 of the FSM: stopping state for the variator.
    case Stop: {
        // variator_m->population()->keepSurvivors(archive_);
        if (dump_offspring_m == false) {
            dumpPopulation(variator_m->population());
        }
        // dump Pareto front
        std::ostringstream filename;
        filename << resultDir_m << "/" << "ParetoFront_" << resultFile_m
                 << "_" << comms_.island_id;
        dumpPopulationToJSON(paretoFront_m, filename, false);

        curState_m = VariatorStopped;

        // notify pilot that we have converged
        int dummy = 0;
        MPI_Request req;
        MPI_Isend(&dummy, 1, MPI_INT, comms_.master_local_pid,
                  MPI_OPT_CONVERGED_TAG, comms_.opt, &req);

        break;
    }

    // State 7 of the FSM: stopping state for the selector.
    // case SelectorStopped: {
    //     curState_m = Stop;
    //     break;
    // }

    // State 8 of the FSM: reset the variator, restart in state 0.
    // case Reset: {
    //     act_gen = 1;
    //     variator_m->population()->keepSurvivors(archive_);
    //     dumpPopulation(variator_m->population());

    //     variator_m->population()->clean_population();
    //     curState_m = ReadyForReset;
    //     break;
    // }

    // State 11 of the FSM: selector has just reset and is ready to
    // start again in state 1.
    // case Restart: {
    //     curState_m = Initialize;
    //     break;
    // }


//-----------------------|    selector STM   |-----------------------------//

    case InitializeSelector: {

        selection();

        // write arc file (all individuals that could ever be used again)
        typename std::map<unsigned int, individual >
            ::iterator it;
        for(it =  variator_m->population()->begin();
            it != variator_m->population()->end(); it++) {
            //archive_.insert(it->first);
            pp_all.push_back(it->first);
        }

        curState_m = Variate;

        break;
    }

    case Select: {
        selection();
        curState_m = Variate;

        break;
    }

    // variator terminated
    case VariatorStopped: {
        curState_m = VariatorTerminate;
        break;
    }

    // // variator ready for reset
    // case ReadyForReset: {
    //     curState_m = ReadyForResetS;
    //     break;
    // }

    // // reset
    // case ReadyForResetS: {
    //     curState_m = Restart;
    //     break;
    // }

    case VariatorTerminate:
    default:
        break;

    }
}

template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::dumpPopulation(boost::shared_ptr<Population_t> population) {
    std::ostringstream filename;
    int fileNumber = act_gen;
    if (dump_offspring_m == false) fileNumber--; // parents are from generation earlier (keeping filenumbers the same)
    filename << resultDir_m << "/" << fileNumber << "_" << resultFile_m
             << "_" << comms_.island_id;

    // only dump old data format if the user requests it
    if (dump_dat_m == true)
        dumpPopulationToFile(population, filename, dump_offspring_m);

    dumpPopulationToJSON(population, filename, dump_offspring_m);
}

template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::dumpPopulationToFile(boost::shared_ptr<Population_t> population,
                                                  std::ostringstream& filename,
                                                  bool dump_offspring) {

    std::ofstream file;
    file.open(filename.str().c_str(), std::ios::out);

    typename std::map<unsigned int, individual >::iterator it;
    it = population->begin();
    // find maximum length of ID
    auto maxID = it->first;
    for(it ++; it != population->end(); it++) {
        if (it->first > maxID)
            maxID = it->first;
    }
    const size_t numDigits = std::to_string(maxID).length();
    size_t next = 0, last = 0;
    std::string delimiter = ",";
    next = file_param_descr_.find(delimiter, last);
    file << std::setw(numDigits + 1) << file_param_descr_.substr(last, next - last) << " ";
    last = next + 1;
    while ((next = file_param_descr_.find(delimiter, last)) != std::string::npos) {
        size_t next2 = file_param_descr_.substr(last, next - last).find(":");
        if (next2 != std::string::npos) {
            last += next2 + 1;
            file << "DVAR: ";
        }
        file << std::setw(14) << file_param_descr_.substr(last, next - last) << " ";
        last = next + 1;
    }
    file << std::setw(14) << file_param_descr_.substr(last) << std::endl;

    file.precision(6);
    file << std::scientific;

    if (dump_offspring == true) {
        for (auto id : finishedBuffer_m) {
            auto ind = population->get_staging(id); // get from staging area
            if (ind) // pointer might be uninitialised (should not happen)
                dumpIndividualToFile(id, ind, file, numDigits);
            else
                std::cout << "Individual " << id << " from buffer not found!" << std::endl;
        }
    } else {
        for(it  = population->begin(); it != population->end(); it++) {
            dumpIndividualToFile(it->first, it->second, file, numDigits);
        }
    }

    file.close();
}

template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::dumpIndividualToFile(int idx,
                                                  individual& ind,
                                                  std::ofstream& file,
                                                  const size_t numDigits) {

    file << std::setw(numDigits + 1) << idx << " ";

    for(size_t i=0; i<ind->objectives_m.size(); i++)
        file << std::setw(14) << ind->objectives_m[i] << " ";
    file << "      ";
    for(size_t i=0; i<ind->genes_m.size(); i++)
        file << std::setw(14) << ind->genes_m[i] << " ";

    file << std::endl;
}

template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::dumpPopulationToJSON(boost::shared_ptr<Population_t> population,
                                                  std::ostringstream& filename,
                                                  bool dump_offspring) {

    typedef boost::property_tree::ptree ptree_t;
    ptree_t tree;

    tree.put("name", "opt-pilot");
    tree.put(OPAL_PROJECT_NAME " version", OPAL_PROJECT_VERSION);
    tree.put("git revision", Util::getGitRevision());

    std::stringstream bounds;
    DVarContainer_t::iterator itr = dvars_m.begin();
    for (bounds_t::iterator it = dVarBounds_m.begin();
         it != dVarBounds_m.end(); ++it, ++itr)
    {
        std::string dvar = boost::get<VAR_NAME>(itr->second);
        bounds << "[ " << it->first << ", " << it->second << " ]";
        tree.put("dvar-bounds." + dvar, bounds.str());
        bounds.str("");
    }

    ptree_t constraints;

    for ( Expressions::Named_t::iterator it = constraints_m.begin();
          it != constraints_m.end(); ++it )
    {
        std::string s = it->second->toString();
        /// cleanup string to make json reader happy
        s.erase(std::remove(s.begin(), s.end(), '"'), s.end());

        // 22. November 2018
        // https://stackoverflow.com/questions/2114466/creating-json-arrays-in-boost-using-property-trees
        ptree_t constraint;
        constraint.put("", s);
        constraints.push_back(std::make_pair("", constraint));
    }

    tree.add_child("constraints", constraints);

    if (dump_offspring == true) {
        for (auto id : finishedBuffer_m) {
            auto ind = population->get_staging(id); // get from staging area
            if (ind) // pointer might be uninitialised (should not happen)
                dumpIndividualToJSON(id, ind, tree);
            else
                std::cout << "Individual " << id << " from buffer not found!" << std::endl;
        }
    } else {
        for (auto it = population->begin(); it != population->end(); it++) {
            dumpIndividualToJSON(it->first, it->second, tree);
        }
    }

    filename << ".json";
    boost::property_tree::write_json(filename.str(), tree);
}

template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::dumpIndividualToJSON(int idx,
                                                  individual& ind,
                                                  boost::property_tree::ptree& tree) {

    std::string id = std::to_string(idx);

    Expressions::Named_t::iterator expr_it;
    expr_it = objectives_m.begin();

    for(size_t i=0; i < ind->objectives_m.size(); i++, expr_it++) {
        std::string name = expr_it->first;
        tree.put("population." + id + ".obj." + name, ind->objectives_m[i]);
    }

    size_t i = 0;
    for(DVarContainer_t::iterator itr = dvars_m.begin(); itr != dvars_m.end(); ++i, ++itr) {
        std::string name = boost::get<VAR_NAME>(itr->second);
        tree.put("population." + id + ".dvar." + name, ind->genes_m[i]);
    }
}


/*-----------------------| selection functions|--------------------------*/

template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::selection() {

    /* Join offspring individuals from variator to population */
    mergeOffspring();

    int size = pp_all.size();

    /* Create internal data structures for selection process */
    copies.resize(size);
    dist.resize(size);
    for (int i = 0; i < size; i++) front.push_back(std::vector<int>(size));


    /* Calculates NSGA2 fitness values for all individuals */
    calcFitnesses();

    /* Calculates distance cuboids */
    calcDistances();

    /* Performs environmental selection
       (truncates 'pp_all' to size 'alpha') */
    environmentalSelection();

    /* Performs mating selection
       (fills mating pool / offspring population pp_sel */
    matingSelection();


    copies.clear();
    dist.clear();
    front.clear();
}


template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::mergeOffspring() {
    pp_all.insert(pp_all.end(), to_selector_.begin(), to_selector_.end());
    to_selector_.clear();
}


template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::calcFitnesses()
{
    int i, j, l;
    int size = pp_all.size();

    std::vector<int> d(size,0); // 1 if dominated by another ind, 0 if not
    std::vector<int> f(size,1); // 1 if not yet in a front, -1 if in a front

    /* initialize fitness and strength values */
    for (i = 0; i < size; i++) {
        fitness_.insert(std::pair<size_t, double>(pp_all[i], 0.0));
        copies[i] = 0;
    }

    /* calculate strength values */
    int num = size; // number of individuals not yet in a front
    for (l = 0; l < size; l++) { // loop over fronts, there can be maximally `size` fronts
        /* find next front */
        for (i = 0; i < size; i++) {
            d[i] = 0; // reset dominate flag
            if (f[i] == -1) continue; // don't consider if already in a front
            for (j = 0; j < size && j != i; j++) {
                if (f[j] == -1) continue;
                individual ind_i = variator_m->population()->get_individual(pp_all[i]);
                individual ind_j = variator_m->population()->get_individual(pp_all[j]);
                if (dominates(ind_j, ind_i)) {
                    d[i] = 1;
                    break;
                }
            }
        }

        /* extract front */
        for (i = 0; i < size; i++) {
            if (f[i] != -1 && d[i] == 0) {
                // add to front
                fitness_[pp_all[i]] = l;
                f[i] = -1;
                num--;
                front[l][copies[l]] = i;
                copies[l] += 1;
            }
        }

        if (num == 0)
            break;
    }
}


template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::calcDistances()
{
    int i, j, l;
    int size = pp_all.size();
    double dmax = 1E99 / (dim_m + 1);

    for (i = 0; i < size; i++) {
        dist[i] = 1;
    }

    for (l = 0; l < size; l++) {
        for (size_t d = 0; d < dim_m; d++) {
            /* sort accorting to d-th objective */
            for (i = 0; i < copies[l]; i++) {
                int min_index = -1;
                int min = i;
                size_t idx1 = pp_all[front[l][i]];
                individual ind1 =
                    variator_m->population()->get_individual(idx1);
                double obj_min = ind1->objectives_m[d];

                for (j = i + 1; j < copies[l]; j++) {

                    // size_t idx1 = pp_all[front[l][j]];
                    size_t idx2 = pp_all[front[l][j]];
                    // individual ind1 =
                    //     variator_m->population()->get_individual(idx1);
                    individual ind2 =
                        variator_m->population()->get_individual(idx2);
                    // double obj1 = ind1->objectives_m[d];
                    double obj2 = ind2->objectives_m[d];

                    if ( obj2 < obj_min ) {
                        min = j;
                        obj_min = obj2;
                    }
                }
                min_index = front[l][min];
                front[l][min] = front[l][i];
                front[l][i] = min_index;
            }

            /* add distances */
            for (i = 0; i < copies[l]; i++) {
                if (i == 0 || i == copies[l] - 1)
                    dist[front[l][i]] += dmax;
                else {
                    size_t idx1 = pp_all[front[l][i+1]];
                    size_t idx2 = pp_all[front[l][i-1]];
                    individual ind1 =
                        variator_m->population()->get_individual(idx1);
                    individual ind2 =
                        variator_m->population()->get_individual(idx2);
                    double obj1 = ind1->objectives_m[d];
                    double obj2 = ind2->objectives_m[d];

                    dist[front[l][i]] += obj1 - obj2;
                }
            }
        }
    }
}


template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::environmentalSelection() {

    int size = pp_all.size();

    for (int i = 0; i < size; i++)
        fitness_[pp_all[i]] += 1.0 / dist[i];

    // get alpha_m fittest individuals
    for (size_t i = 0; i < alpha_m; i++) {
        int min = i;
        for (int j = i + 1; j < size; j++) {
            if (fitness_[pp_all[j]] < fitness_[pp_all[min]])
                min = j;
        }
        // swap
        std::swap(pp_all[min], pp_all[i]);
    }
    // erase others
    pp_all.erase(pp_all.begin() + alpha_m, pp_all.end());
}


/* Fills mating pool 'parent_queue_' */
template< template <class> class CO, template <class> class MO >
void FixedPisaNsga2<CO, MO>::matingSelection() {
    // select lambda_m parents
    for (size_t i = 0; i < lambda_m; i++) {
        int winner = irand(pp_all.size());
        for (int j = 0; j < tournament_m; j++) {
            int opponent = irand(pp_all.size());
            if (fitness_[pp_all[opponent]] < fitness_[pp_all[winner]]
                || winner == opponent) {
                winner = opponent;
            }
        }
        parent_queue_.push_back(pp_all[winner]);
    }
}


/* Determines if one individual dominates another.
   Minimizing fitness values. */
template< template <class> class CO, template <class> class MO >
int FixedPisaNsga2<CO, MO>::dominates(individual ind_a, individual ind_b) {

    int a_is_worse = 0, b_is_worse = 0;
    // int equal = 1;

    for (size_t i = 0; i < ind_a->objectives_m.size()/* && !a_is_worse*/; i++) {
        if (ind_a->objectives_m[i] > ind_b->objectives_m[i]) a_is_worse = 1;
        if (ind_a->objectives_m[i] < ind_b->objectives_m[i]) b_is_worse = 1;
        // a_is_worse = ind_a->objectives_m[i] > ind_b->objectives_m[i];
        // equal = (ind_a->objectives_m[i] == ind_b->objectives_m[i]) && equal;
    }

    return (b_is_worse && !a_is_worse);
    // return (!equal && !a_is_worse);
}

template< template <class> class CO, template <class> class MO >
bool FixedPisaNsga2<CO, MO>::checkParetoFront(unsigned int newid) {
    individual newInd = variator_m->population()->get_staging(newid);

    for (auto it = paretoFront_m->begin(); it != paretoFront_m->end(); it++) {
        if (dominates(it->second, newInd) == true) return false;
    }
    //check if same (should happen rarely)
    if (paretoFront_m->isRepresentedInPopulation(newInd->genes_m) == true)
        return false;

    for (auto it = paretoFront_m->begin(); it != paretoFront_m->end(); it++) {
        if (dominates(newInd, it->second) == true) {
            // remove individuals no longer in Pareto front
            it = paretoFront_m->erase(it);
        }
    }
    paretoFront_m->add_individual(newInd);
    paretoFront_m->commit_individuals();

    return true;
}

template< template <class> class CO, template <class> class MO >
std::string FixedPisaNsga2<CO, MO>::getStateString(PisaState_t state) const {
    switch(state) {
    case Initialize:
        return "Initialize";
    case InitializeSelector:
        return "InitializeSelector";
    case Variate:
        return "Variate";
    case Select:
        return "Select";
    case Stop:
        return "Stop";
    case VariatorStopped:
        return "VariatorStopped";
    case VariatorTerminate:
        return "VariatorTerminate";
    // case SelectorStopped:
    //     return "SelectorStopped";
    // case Reset:
    //     return "Reset";
    // case ReadyForReset:
    //     return "ReadyForReset";
    // case ReadyForResetS:
    //     return "ReadyForResetS";
    // case Restart:
    //     return "Restart";
    default:
        return "";
    }
}