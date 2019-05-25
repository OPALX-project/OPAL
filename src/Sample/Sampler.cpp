#include <iostream>
#include <string>
#include <limits>

#include "Sample/Sampler.h"

#include "OPALconfig.h"
#include "Utilities/Util.h"

#include "Util/OptPilotException.h"
#include "Util/MPIHelper.h"

#include <boost/property_tree/json_parser.hpp>

#include <boost/filesystem.hpp>

Sampler::Sampler(Expressions::Named_t objectives,
                 Expressions::Named_t constraints,
                 DVarContainer_t dvars,
                 size_t dim, Comm::Bundle_t comms,
                 CmdArguments_t args,
                 std::vector<double> hypervolRef,
                 int nrWorkerGroups)
    : Optimizer(comms.opt)
{
    throw OptPilotException("Sampler::Sampler",
                            "We shouldn't get here!");
}


Sampler::Sampler(const std::map<std::string,
                                std::shared_ptr<SamplingMethod>
                    >& sampleMethods,
                 Expressions::Named_t objectives,
                 DVarContainer_t dvars,
                 Comm::Bundle_t comms,
                 CmdArguments_t args)
    : Optimizer(comms.opt)
    , sampleMethods_m(sampleMethods)
    , comms_(comms)
    , dvars_m(dvars)
    , objectives_m(objectives)
    , args_(args)
{
    my_local_pid_ = 0;
    MPI_Comm_rank(comms_.opt, &my_local_pid_);
    
    std::string fname = "samples_" + std::to_string(comms_.island_id) + ".json";

    resultFile_m = args->getArg<std::string>("outfile", fname, false);
    resultDir_m = args->getArg<std::string>("outdir", "samples", false);

    if ( !boost::filesystem::exists(resultDir_m) ) {
        boost::filesystem::create_directory(resultDir_m);
    }

    DVarContainer_t::iterator itr;
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        dVarBounds_m.push_back(
                std::pair<double, double>
                    (boost::get<LOWER_BOUND>(itr->second),
                     boost::get<UPPER_BOUND>(itr->second)));
    }
}


void Sampler::initialize() {

    nSamples_m = args_->getArg<int>("nsamples", true);
    act_sample_m = 0;
    done_sample_m = 0;
    curState_m = SUBMIT;

    int nMasters = args_->getArg<int>("num-masters", true);

    if ( nMasters > nSamples_m )
        throw OptPilotException("Sampler::initialize",
                                "More masters than samples.");
    
    // unique job id
    int nLocSamples = nSamples_m / nMasters;
    int rest = nSamples_m - nMasters * nLocSamples;
    
    if ( comms_.island_id < rest )
        nLocSamples++;
    
    if ( rest == 0 )
        gid = nLocSamples * comms_.island_id;
    else {
        if ( comms_.island_id < rest ) {
            gid = nLocSamples * comms_.island_id;
        } else {
            gid = (nLocSamples + 1) * rest + (comms_.island_id - rest) * nLocSamples;
        }
    }
    
    nSamples_m = nLocSamples;
    
    // start poll loop
    run();
}


bool Sampler::onMessage(MPI_Status status, size_t length) {
    MPITag_t tag = MPITag_t(status.MPI_TAG);
    switch(tag) {
        case REQUEST_FINISHED: {
            unsigned int jid = static_cast<unsigned int>(length);
            typename std::map<size_t, boost::shared_ptr<Individual_t> >::iterator it;
            it = jobmapping_m.find(jid);

            if(it == jobmapping_m.end()) {
                std::cout << "NON-EXISTING JOB with ID = " << jid << std::endl;
                throw OptPilotException("Sampler::onMessage",
                        "non-existing job");
            }


            boost::shared_ptr<Individual_t> ind = it->second;

            reqVarContainer_t res;
            MPI_Recv_reqvars(res, status.MPI_SOURCE, comms_.opt);

            ind->objectives.clear();

            reqVarContainer_t::iterator itr = res.begin();
            for(; itr != res.end(); ++itr) {
                // mark invalid if expression could not be evaluated or constraint does not hold
                if(!itr->second.is_valid || (itr->second.value.size() > 1 && !itr->second.value[0])) {
                    ind->objectives.push_back(std::numeric_limits<double>::infinity());
                } else {
                    // update objective value for valid objective
                    if(itr->second.value.size() == 1)
                        ind->objectives.push_back(itr->second.value[0]);
                }
            }

            addIndividualToJSON_m(ind);

            jobmapping_m.erase(it);

            done_sample_m++;

            return true;
        }
        default: {
            std::cout << "(Sampler) Error: unexpected MPI_TAG: "
                      << status.MPI_TAG << std::endl;
            return false;
        }
    }
}


void Sampler::postPoll() {

    if ( act_sample_m < nSamples_m ) {
        this->createNewIndividual_m();
    }

    runStateMachine();
}


void Sampler::createNewIndividual_m() {

    std::vector<std::string> dNames;

    DVarContainer_t::iterator itr;
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        std::string dName = boost::get<VAR_NAME>(itr->second);
        dNames.push_back(dName);
    }

    boost::shared_ptr<Individual_t> ind = boost::shared_ptr<Individual_t>( new Individual_t(dNames));

    ind->id = gid++;
    
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        std::string dName = boost::get<VAR_NAME>(itr->second);
        int i = ind->getIndex(dName);
        sampleMethods_m[dName]->create(ind, i);
    }


    individuals_m.push(ind);
}

void Sampler::dumpIndividualsToJSON_m() {

    tree_m.put("name", "sampler");
    tree_m.put(OPAL_PROJECT_NAME " version", OPAL_PROJECT_VERSION);
    tree_m.put("git revision", Util::getGitRevision());

    std::stringstream bounds;
    DVarContainer_t::iterator itr = dvars_m.begin();
    for (bounds_t::iterator it = dVarBounds_m.begin();
         it != dVarBounds_m.end(); ++it, ++itr)
    {
        std::string dvar = boost::get<VAR_NAME>(itr->second);
        bounds << "[ " << it->first << ", " << it->second << " ]";
        tree_m.put("dvar-bounds." + dvar, bounds.str());
        bounds.str("");
    }

    std::ostringstream filename;
    filename << resultDir_m << "/" << resultFile_m
             << "_samples_" << comms_.island_id << ".json";

    boost::property_tree::write_json(filename.str(), tree_m);
}

void Sampler::addIndividualToJSON_m(const boost::shared_ptr<Individual_t>& ind) {
    
    std::string id = std::to_string(ind->id);
    
    DVarContainer_t::iterator itr;
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        std::string name = boost::get<VAR_NAME>(itr->second);
        int i = ind->getIndex(name);
        tree_m.put("samples." + id + ".dvar." + name, ind->genes[i]);
    }

    Expressions::Named_t::iterator expr_it;
    expr_it = objectives_m.begin();

    for(size_t i=0; i < ind->objectives.size(); i++, expr_it++) {
        std::string name = expr_it->first;

        // skip dummy objective (SamplePilot.h, line 64)
        if ( name == "dummy" )
            continue;

        tree_m.put("samples." + id + ".obj." + name, ind->objectives[i]);
    }
}


void Sampler::runStateMachine() {

    switch(curState_m) {

        case SUBMIT: {
            if ( done_sample_m == nSamples_m) {
                curState_m = STOP;
            } else {
                if ( act_sample_m != nSamples_m ) {
                    dispatch_forward_solves();
                }
            }
            break;
        }
        case STOP: {

            dumpIndividualsToJSON_m();

            curState_m = TERMINATE;

            // notify pilot that we have converged
            int dummy = 0;
            MPI_Request req;
            MPI_Isend(&dummy, 1, MPI_INT, comms_.master_local_pid,
                      MPI_OPT_CONVERGED_TAG, comms_.opt, &req);

            break;
        }

        case TERMINATE: {
            break;
        }
    }
}


void Sampler::dispatch_forward_solves() {

    while ( !individuals_m.empty() ) {
        boost::shared_ptr<Individual_t> ind = individuals_m.front();

        individuals_m.pop();

        Param_t params;
        DVarContainer_t::iterator itr;

        for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
            std::string dName = boost::get<VAR_NAME>(itr->second);
            int i = ind->getIndex(dName);
            params.insert(
                std::pair<std::string, double>
                    (dName, ind->genes[i]));
        }

        size_t jid = static_cast<size_t>(ind->id);
        int pilot_rank = comms_.master_local_pid;


        act_sample_m++;

        // now send the request to the pilot
        MPI_Send(&jid, 1, MPI_UNSIGNED_LONG, pilot_rank, OPT_NEW_JOB_TAG, comms_.opt);

        MPI_Send_params(params, pilot_rank, comms_.opt);

        jobmapping_m.insert(
                std::pair<size_t, boost::shared_ptr<Individual_t> >(jid, ind));
    }
}