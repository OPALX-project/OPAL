#include <iostream>
#include <sstream>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cstdlib>
#include <vector>
#include <ctime>

#include "Simulation/OpalSimulation.h"

#include "Util/SDDSReader.h"
#include "Util/OptPilotException.h"
#include "Util/HashNameGenerator.h"

#include "Expression/SumErrSq.h"
#include "Expression/FromFile.h"

#include "boost/variant.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/algorithm/string.hpp"

#include "boost/filesystem.hpp"
#include "boost/filesystem/operations.hpp"

// access to OPAL lib
#include "opal.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

OpalSimulation::OpalSimulation(Expressions::Named_t objectives,
                               Expressions::Named_t constraints,
                               Param_t params, std::string name,
                               MPI_Comm comm, int groupId,
                               CmdArguments_t args)
               : Simulation(args)
               , objectives_(objectives)
               , constraints_(constraints)
               , comm_(comm)
{
    namespace fs = boost::filesystem;

    if(getenv("SIMTMPDIR") == NULL) {
        std::cout << "Environment variable SIMTMPDIR not defined!"
                  << std::endl;
        simTmpDir_ = getenv("PWD");
    } else
        simTmpDir_ = getenv("SIMTMPDIR");

    simulationName_ = name;

    // prepare design variables given by the optimizer for generating the
    // input file
    std::vector<std::string> dict;
    for(auto parameter : params) {
        std::ostringstream tmp;
        tmp.precision(15);
        tmp << parameter.first << "=" << parameter.second;
        dict.push_back(tmp.str());

        std::ostringstream value;
        value.precision(15);
        value << parameter.second;
        userVariables_.insert(
            std::pair<std::string, std::string>(parameter.first, value.str()));
    }

    // hash the dictionary to get a short unique directory name for temporary
    // simulation data
    std::string hash = HashNameGenerator::generate(dict);
    std::ostringstream tmp;
    tmp.precision(15);

    tmp << simTmpDir_ << "/" << hash << "_" << groupId;

    simulationDirName_ = tmp.str();

    if(getenv("TEMPLATES") == NULL) {
        throw OptPilotException("OpalSimulation::OpalSimulation",
            "Environment variable TEMPLATES not defined!");
    }
    std::string tmplDir = getenv("TEMPLATES");
    std::string tmplFile = tmplDir + "/" + simulationName_ + ".tmpl";
    // data file is assumed to be located in the root directory
    std::string dataFile = simulationName_ + ".data";
    if (!fs::exists(dataFile))
        throw OptPilotException("OpalSimulation::OpalSimulation", "The data file '" + dataFile + "' doesn't exist");

    if (!fs::exists(tmplFile))
        throw OptPilotException("OpalSimulation::OpalSimulation", "The template file '" + tmplFile + "' doesn't exit");

    gs_.reset(new GenerateOpalSimulation(tmplFile, dataFile, userVariables_));
}


OpalSimulation::~OpalSimulation() {
    requestedVars_.clear();
    userVariables_.clear();
}

bool OpalSimulation::hasResultsAvailable() {

    std::string infile = simulationDirName_ + "/" + simulationName_ + ".in";
    struct stat fileInfo;

    if(stat(infile.c_str(), &fileInfo) == 0) {
        std::cout << "-> Simulation input file (" << infile
                  << ") already exist from previous run.." << std::endl;
        return true;
    }

    return false;
}


void OpalSimulation::setupSimulation() {

    // only on processor in comm group has to setup files
    int rank = 0;
    MPI_Comm_rank(comm_, &rank);
    if(rank == 0) {

        mkdir((const char*)(simulationDirName_.c_str()), 0777);

        std::string infile = simulationDirName_ + "/" +
                             simulationName_ + ".in";
        gs_->writeInputFile(infile);

        // linking fieldmaps
        if(getenv("FIELDMAPS") == NULL) {
            throw OptPilotException("OpalSimulation::OpalSimulation",
                "Environment variable FIELDMAPS not defined!");
        }
        std::string fieldmapPath = getenv("FIELDMAPS");

        struct dirent **files;
        int count = scandir(fieldmapPath.c_str(), &files, 0, alphasort);

        for(int i=0; i<count; i++) {
            if (files[i]->d_name == std::string(".") ||
                files[i]->d_name == std::string("..")) continue;
            std::string source = fieldmapPath + "/" + files[i]->d_name;
            std::string target = simulationDirName_ + '/' + files[i]->d_name;
	    int err = symlink(source.c_str(), target.c_str());
	    if (err != 0) {
	      // FIXME properly handle error
	      std::cout << "Cannot symlink fieldmap "
			<< source.c_str() << " to "
			<< target.c_str() << " error no " << err << std::endl;
	      std::cout << "fieldmapPath " << fieldmapPath << " i= " << i << std::endl;
	      std::cout << "target       " << simulationDirName_ + '/' + files[i]->d_name << std::endl;

	    }
	}
    }

    MPI_Barrier(comm_);
}


void OpalSimulation::redirectOutToFile() {

    // backup stdout and err file handles
    strm_buffer_ = std::cout.rdbuf();
    strm_err_ = std::cerr.rdbuf();

    int world_pid = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_pid);

    std::ostringstream fname;
    fname << "sim.out." << world_pid;
    std::ofstream file(fname.str().c_str());
    fname << ".err";
    std::ofstream err(fname.str().c_str());

    // and redirect stdout and err to new files
    std::cout.rdbuf(file.rdbuf());
    std::cerr.rdbuf(err.rdbuf());
}


void OpalSimulation::restoreOut() {
    std::cout.rdbuf(strm_buffer_);
    std::cerr.rdbuf(strm_err_);
}


void OpalSimulation::run() {

    // make sure input file is not already existing
    MPI_Barrier(comm_);
    if( hasResultsAvailable() ) return;
    MPI_Barrier(comm_);

    setupSimulation();

    pwd_ = getenv("PWD");
    pwd_ += "/";
    int err = chdir(simulationDirName_.c_str());

    if (err != 0) {
        std::cout << "Cannot chdir to "
                  << simulationDirName_.c_str() << std::endl;
        std::cout << "Continuing 1, disregarding this simulation.."
                  << std::endl;
        return;
    }

    // setup OPAL command line options
    std::ostringstream inputFileName;
    inputFileName << simulationName_ << ".in";

    char exe_name[] = "opal";
    char *inputfile = new char[inputFileName.str().size()+1] ;
    strcpy(inputfile, inputFileName.str().c_str());
    char nocomm[] = "--nocomminit";
    char info[] = "--info";
    char info0[] = "0";
    char warn[] = "--warn";
    char warn0[] = "0";
    char *arg[] = { exe_name, inputfile, nocomm, info, info0, warn, warn0 };

    int seed = Options::seed;

    try {

        //FIXME: this seems to crash OPAL in some cases
        //redirectOutToFile();
#ifdef SUPRESS_OUTPUT
        //XXX: hack to disable output to stdout and stderr
        std::cout.setstate(std::ios::failbit);
        // std::cerr.setstate(std::ios::failbit);
#endif

        // now we can run the simulation
        run_opal(arg, inputFileName.str(), -1, comm_);

        //restoreOut();
#ifdef SUPRESS_OUTPUT
        std::cout.clear();
        std::cerr.clear();
#endif

    } catch(OpalException *ex) {

        //restoreOut();
#ifdef SUPRESS_OUTPUT
        std::cout.clear();
        std::cerr.clear();
#endif

        std::cout << "Opal exception during simulation run: "
                  << ex->what() << std::endl;
        std::cout << "Continuing 2, disregarding this simulation.."
                  << std::endl;

    } catch(ClassicException *ex) {

        //restoreOut();
#ifdef SUPRESS_OUTPUT
        std::cout.clear();
        std::cerr.clear();
#endif

        std::cout << "Classic exception during simulation run: "
                  << ex->what() << std::endl;
        std::cout << "Continuing 3, disregarding this simulation.."
                  << std::endl;
    }

    Options::seed = seed;

    delete[] inputfile;
    err = chdir(pwd_.c_str());
    if (err != 0) {
        std::cout << "Cannot chdir to "
                  << simulationDirName_.c_str() << std::endl;
    }
}


void OpalSimulation::collectResults() {

    // clear old solutions
    requestedVars_.clear();

    int err = chdir(simulationDirName_.c_str());
    if (err != 0) {
        std::cout << "Cannot chdir to "
                  << simulationDirName_.c_str() << std::endl;
        std::cout << "Continuing, with cleanup.."
                  << std::endl;
        cleanUp();
        return;
    }

    std::string fn = simulationName_ + ".stat";
    struct stat fileInfo;

    // if no stat file, simulation parameters produced invalid bunch
    if(stat(fn.c_str(), &fileInfo) != 0) {
        invalidBunch();
    } else {

        for(auto namedObjective : objectives_) {

            Expressions::Expr_t *objective = namedObjective.second;

            // find out which variables we need in order to evaluate the
            // objective
            variableDictionary_t variable_dictionary;
            bool check = getVariableDictionary(variable_dictionary,fn,objective);
            if (check == false) break;

            // and evaluate the expression using the built dictionary of
            // variable values
            Expressions::Result_t result =
                objective->evaluate(variable_dictionary);

            std::vector<double> values;
            values.push_back(boost::get<0>(result));
            bool is_valid = boost::get<1>(result);

            reqVarInfo_t tmps = {EVALUATE, values, is_valid};
            requestedVars_.insert(
                std::pair<std::string, reqVarInfo_t>(namedObjective.first, tmps));

        }

        // .. and constraints
        for(auto namedConstraint : constraints_) {

            Expressions::Expr_t *constraint = namedConstraint.second;

            // find out which variables we need in order to evaluate the
            // objective
            variableDictionary_t variable_dictionary;
            bool check = getVariableDictionary(variable_dictionary,fn,constraint);
            if (check == false) break;

            Expressions::Result_t result =
                constraint->evaluate(variable_dictionary);

            std::vector<double> values;
            values.push_back(boost::get<0>(result));
            bool is_valid = boost::get<1>(result);

            //FIXME: hack to give feedback about values of LHS and RHS
            std::string constr_str = constraint->toString();
            std::vector<std::string> split;
            boost::split(split, constr_str, boost::is_any_of("<>!="),
                        boost::token_compress_on);
            std::string lhs_constr_str = split[0];
            std::string rhs_constr_str = split[1];
            boost::trim_left_if(rhs_constr_str, boost::is_any_of("="));

            functionDictionary_t funcs = constraint->getRegFuncs();
            boost::scoped_ptr<Expressions::Expr_t> lhs(
                new Expressions::Expr_t(lhs_constr_str, funcs));
            boost::scoped_ptr<Expressions::Expr_t> rhs(
                new Expressions::Expr_t(rhs_constr_str, funcs));

            Expressions::Result_t lhs_res = lhs->evaluate(variable_dictionary);
            Expressions::Result_t rhs_res = rhs->evaluate(variable_dictionary);

            values.push_back(boost::get<0>(lhs_res));
            values.push_back(boost::get<0>(rhs_res));

            reqVarInfo_t tmps = {EVALUATE, values, is_valid};
            requestedVars_.insert(
                    std::pair<std::string, reqVarInfo_t>(namedConstraint.first, tmps));

        }

    }

    err = chdir(pwd_.c_str());
    if (err != 0) {
        std::cout << "Cannot chdir to "
                  << simulationDirName_.c_str() << std::endl;
    }

    cleanUp();
}

bool OpalSimulation::getVariableDictionary(variableDictionary_t& dictionary,
                                           const std::string& filename,
                                           const Expressions::Expr_t* const expression) {

    std::set<std::string> req_vars = expression->getReqVars();
    if(req_vars.empty()) return true;

    boost::scoped_ptr<SDDSReader> sddsr(new SDDSReader(filename));
    try {
        sddsr->parseFile();
    } catch(OptPilotException &e) {
        std::cout << "Exception while parsing SDDS file: "
                  << e.what() << std::endl;

        //XXX: in this case we mark the bunch as invalid since
        //     broken stat files can crash opt (why do they
        //     exist?)
        invalidBunch();
        return false;
    }

    // get all the required variable values from the stat file
    for(std::string req_var : req_vars) {
        if(dictionary.count(req_var) != 0) continue;

        try {
            double value = 0.0;
            sddsr->getValue(-1 /*atTime*/, req_var, value);
            dictionary.insert(std::pair<std::string, double>(req_var, value));
        } catch(OptPilotException &e) {
            std::cout << "Exception while getting value "
                      << "from SDDS file: " << e.what()
                      << std::endl;
            return false;
        }
    }
    return true;
}

void OpalSimulation::invalidBunch() {

    for(auto namedObjective : objectives_) {
        std::vector<double> tmp_values;
        tmp_values.push_back(0.0);
        reqVarInfo_t tmps = {EVALUATE, tmp_values, false};
        requestedVars_.insert(
                std::pair<std::string, reqVarInfo_t>(namedObjective.first, tmps));
    }
}

void OpalSimulation::cleanUp() {

#ifdef BOOST_FILESYSTEM
    try {
        boost::filesystem::path p(simulationDirName_.c_str());
        boost::filesystem::remove_all(p);
    } catch(...) {
    }
#endif

}
