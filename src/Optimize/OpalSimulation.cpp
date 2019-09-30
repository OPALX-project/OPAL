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
#include <exception>

#include "Optimize/OpalSimulation.h"

#include "Util/SDDSReader.h"
#include "Util/SDDSParser.h"
#include "Util/SDDSParser/SDDSParserException.h"
#include "Util/OptPilotException.h"
#include "Util/NativeHashGenerator.h"

#include "AbstractObjects/OpalData.h"

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
#include "Utilities/Util.h"

OpalSimulation::OpalSimulation(Expressions::Named_t objectives,
                               Expressions::Named_t constraints,
                               Param_t params, std::string name,
                               MPI_Comm comm, CmdArguments_t args)
               : Simulation(args)
               , objectives_(objectives)
               , constraints_(constraints)
               , comm_(comm)
               , id_m(-1)
{
    namespace fs = boost::filesystem;

    simTmpDir_ = args->getArg<std::string>("simtmpdir");
    if (simTmpDir_ == "") {
        if(getenv("SIMTMPDIR") == NULL) {
            std::cout << "Environment variable SIMTMPDIR not defined!"
                      << std::endl;
            simTmpDir_ = getenv("PWD");
        } else
            simTmpDir_ = getenv("SIMTMPDIR");
    }
    simulationName_ = name;

    // prepare design variables given by the optimizer for generating the
    // input file
    std::vector<std::string> dict;
    for(auto parameter : params) {
        std::ostringstream tmp;
        tmp.precision(15);
        tmp << parameter.first << "=" << parameter.second;
        dvarNames_.insert(parameter.first);
        dict.push_back(tmp.str());

        std::ostringstream value;
        value.precision(15);
        value << parameter.second;
        userVariables_.insert(
            std::pair<std::string, std::string>(parameter.first, value.str()));
    }

    /*
      This is a copy from Comm/Splitter/ManyMasterSplit.h
      in order to calculate the leader which is the unique ID in case
      of more than one core per worker.
    */

    int my_rank=0;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int world_size=0;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    unsigned num_coworkers_worker_ = 0;
    num_coworkers_worker_ = args->getArg<size_t>("num-coworkers");

    unsigned group_start = 0;

    unsigned worker_group = ((my_rank % world_size) - 2) / num_coworkers_worker_;

    unsigned leader_ = group_start + 2 + worker_group * num_coworkers_worker_;
    leader_ = leader_ % world_size;

    // hash the dictionary to get a short unique directory name for temporary
    // simulation data
    std::string hash = NativeHashGenerator::generate(dict);

    std::ostringstream tmp;
    tmp.precision(15);

    tmp << hash << "_" << leader_;
    simdir_ = tmp.str();
    tmp.str(std::string());
    tmp << simTmpDir_ << "/" << simdir_;


    simulationDirName_ = tmp.str();

    int nTrials = 100;
    while (fs::exists(simulationDirName_) && nTrials > 0) {
        std::ostringstream tmp1;
        tmp1.precision(15);
        std::string hash1 = NativeHashGenerator::generate(dict);
        tmp1 << hash1 << "_" << leader_;
        simdir_ = tmp1.str();
        tmp1.str(std::string());
        tmp1 << simTmpDir_ << "/" << simdir_;
        simulationDirName_ = tmp1.str();
        --nTrials;
    }

    if (fs:exists(simulationDirName_)) {
        throw OptPilotException("OpalSimulation::OpalSimulation",
                                "Tried to generate a unique directory a 100 times. I give up now.");
    }



    std::string tmplDir = args->getArg<std::string>("templates");
    if (tmplDir == "") {
        if(getenv("TEMPLATES") == NULL) {
            throw OptPilotException("OpalSimulation::OpalSimulation",
                                    "Environment variable TEMPLATES not defined!");
        }
        tmplDir = getenv("TEMPLATES");
    }
    std::string tmplFile = tmplDir + "/" + simulationName_ + ".tmpl";
    // data file is assumed to be located in the root directory
    std::string dataFile = simulationName_ + ".data";
    fs::path pwd = fs::current_path();
    if (!fs::exists(dataFile))
        throw OptPilotException("OpalSimulation::OpalSimulation",
                                "The data file '" + dataFile + "' \n" +
                                "     doesn't exist in directory '" + pwd.native() + "'");

    if (!fs::exists(tmplFile))
        throw OptPilotException("OpalSimulation::OpalSimulation",
                                "The template file '" + tmplFile + "' doesn't exit");

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


void OpalSimulation::createSymlink_m(const std::string& path) {
    namespace fs = boost::filesystem;

    for (auto &p: fs::directory_iterator(path)) {
        fs::path source = p.path();
        fs::path target(simulationDirName_ + "/");
        target +=source.filename();

        try {
            fs::create_symlink(source, target);
        } catch (fs::filesystem_error &e) {
            std::cerr << e.what() << "\n"
                      << "in OpalSimulation::createSymlink()" << std::endl;
        }
    }
}


void OpalSimulation::copyH5_m() {
    CmdArguments_t args = getArgs();
    std::string restartfile = args->getArg<std::string>("restartfile", "", false);

    if (restartfile.empty()) return;

    namespace fs = boost::filesystem;
    if ( !fs::exists(restartfile) ) {
        std::cerr << "H5 file '" + restartfile + "' doesn't exist." << "\n"
                  << "in OpalSimulation::copyH5_m()" << std::endl;

        return;
    }

    try {
        fs::path srcfile(restartfile);
        fs::path targetfile(simulationDirName_ + "/" + simulationName_ + ".h5");
        fs::copy_file(srcfile, targetfile);
    } catch (fs::filesystem_error &ex) {
        std::cerr << ex.what() << "\n"
                  << "in OpalSimulation::copyH5_m()" << std::endl;
    }
}


void OpalSimulation::setupSimulation() {
    namespace fs = boost::filesystem;

    CmdArguments_t args = getArgs();
    std::string restartfile = args->getArg<std::string>("restartfile", "", false);

    if ( id_m > -1 ) {
        std::ostringstream tmp;
        tmp << simTmpDir_ << "/" << id_m;
        simulationDirName_ = tmp.str();
    }
    std::string dataDir = simulationDirName_ + "/data";

    OpalData *opal = OpalData::getInstance();
    opal->setOptimizerFlag();

    // linking fieldmaps + distributions
    if (getenv("FIELDMAPS") == NULL) {
        throw OptPilotException("OpalSimulation::setupSimulation",
                                "Environment variable FIELDMAPS not defined!");
    }

    setupFSStructure();

    MPI_Barrier(comm_);

    if (!fs::exists(simulationDirName_)) {
        throw OptPilotException("OpalSimulation::setupSimulation",
                                "Directory '" + simulationDirName_ + "' doesn't exist");
    }

    if (!fs::exists(dataDir)) {
        throw OptPilotException("OpalSimulation::setupSimulation",
                                "Directory '" + dataDir + "' doesn't exist");
    }

    if (!restartfile.empty() &&
        !fs::exists(simulationDirName_ + "/" + simulationName_ + ".h5")) {
        throw OptPilotException("OpalSimulation::setupSimulation",
                                "H5 file '" + simulationDirName_ + "/" + simulationName_ + ".h5' doesn't exist");
    }
}

void OpalSimulation::setupFSStructure() {
    namespace fs = boost::filesystem;

    int rank = 0;
    MPI_Comm_rank(comm_, &rank);
    if (rank != 0) return;     // only one processor in comm group has to setup files

//     if (fs::exists(simulationDirName_)) {
//         fs::remove_all(simulationDirName_);
//     }

    try {
        fs::create_directory(simulationDirName_);
        fs::permissions(simulationDirName_,
                        fs::owner_all |
                        fs::group_read |
                        fs::group_exe |
                        fs::others_read |
                        fs::others_exe);

    } catch (fs::filesystem_error &e) {
        std::cerr << e.what() << "\n"
                  << "in OpalSimulation::setupSimulation" << std::endl;
        return;
    }

    try {
        std::string dataDir = simulationDirName_ + "/data";

        fs::create_directory(dataDir);
        fs::permissions(dataDir,
                        fs::owner_all |
                        fs::group_read |
                        fs::group_exe |
                        fs::others_read |
                        fs::others_exe);

    } catch (fs::filesystem_error &e) {
        std::cerr << e.what() << "\n"
                  << "in OpalSimulation::setupSimulation" << std::endl;
        return;
    }

    std::string infile = simulationDirName_ + "/" +
        simulationName_ + ".in";
    gs_->writeInputFile(infile);

    std::string fieldmapPath = getenv("FIELDMAPS");
    this->createSymlink_m(fieldmapPath);

    if (getenv("DISTRIBUTIONS") != NULL) {
        std::string distPath = getenv("DISTRIBUTIONS");
        this->createSymlink_m(distPath);
    }

    this->copyH5_m();
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
    namespace fs = boost::filesystem;

    // make sure input file is not already existing
    MPI_Barrier(comm_);
    if( hasResultsAvailable() ) return;
    MPI_Barrier(comm_);

    setupSimulation();

    pwd_ = fs::current_path().native();
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
    char *inputfile = new char[inputFileName.str().size()+1] ;
    strcpy(inputfile, inputFileName.str().c_str());
    int seed = Options::seed;

    CmdArguments_t args = getArgs();
    int restartStep= args->getArg<int>("restartstep",
                                       std::numeric_limits<int>::min(), false);
    std::string restartfile = args->getArg<std::string>("restartfile", "", false);

    try {
        if ( restartStep > -2 && restartfile.empty() ) {
            throw OpalException("OpalSimulation::run()",
                                "Restart specified but no restart H5 file available.");
        }

        char exe_name[] = "opal";
        char nocomm[]   = "--nocomminit";
        char info[]     = "--info";
        char info0[]    = "0";
        char warn[]     = "--warn";
        char warn0[]    = "0";
        char *arg[]     = { exe_name, inputfile, nocomm, info, info0, warn, warn0 };

        //FIXME: this seems to crash OPAL in some cases
        //redirectOutToFile();
#ifdef SUPRESS_OUTPUT
        //XXX: hack to disable output to stdout and stderr
        std::cout.setstate(std::ios::failbit);
        // std::cerr.setstate(std::ios::failbit);
#endif
        // now we can run the simulation
        run_opal(arg, inputFileName.str(), restartStep, Options::infoLevel, Options::warnLevel, comm_);

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

        std::cout << "Opal exception during simulation run: \n"
                  << ex->where() << "\n"
                  << ex->what() << std::endl;
        std::cout << "Continuing 2, disregarding this simulation.."
                  << std::endl;

    } catch(ClassicException *ex) {

        //restoreOut();
#ifdef SUPRESS_OUTPUT
        std::cout.clear();
        std::cerr.clear();
#endif

        std::cout << "Classic exception during simulation run: \n"
                  << ex->where() << "\n"
                  << ex->what() << std::endl;
        std::cout << "Continuing 3, disregarding this simulation.."
                  << std::endl;

    }

    Options::seed = seed;

    delete[] inputfile;
    err = chdir(pwd_.c_str());
    if (err != 0) {
        std::cout << "Cannot chdir to "
                  << pwd_ << std::endl;
    }
}


std::map<std::string, std::vector<double> > OpalSimulation::getData(const std::vector<std::string> &statVariables) {
    std::map<std::string, std::vector<double> > ret;
    SDDS::SDDSParser parser(simulationDirName_ + "/" + simulationName_ + ".stat");
    parser.run();
    for (const std::string &var : statVariables) {
        SDDS::ast::columnData_t column;
        try {
            column = parser.getColumnData(var);
        } catch (SDDSParserException &e) {
            std::cout << "failed to read data: " << e.what() << " in " << e.where() << std::endl;
            continue;
        }

        std::vector<double> values;
        values.reserve(column.size());
        auto type = parser.getColumnType(var);
        for (const auto val: column) {
            values.push_back(parser.getBoostVariantValue<double>(val,(int)type));
        }
        ret.insert(std::make_pair(var, values));
    }

    return ret;
}

void OpalSimulation::collectResults() {

    // std::cout << "collectResults" << std::endl;

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
        Expressions::Named_t::iterator namedIt;
        try {
            for(namedIt=objectives_.begin(); namedIt!=objectives_.end(); ++namedIt) {
                if (namedIt->first == "dummy") continue; // FIXME SamplePilot has default objective named dummy
                Expressions::Expr_t *objective = namedIt->second;

                // find out which variables we need in order to evaluate the objective
                variableDictionary_t variable_dictionary;
                getVariableDictionary(variable_dictionary,fn,objective);

                // and evaluate the expression using the built dictionary of
                // variable values
                Expressions::Result_t result =
                    objective->evaluate(variable_dictionary);

                std::vector<double> values;
                values.push_back(boost::get<0>(result));
                bool is_valid = boost::get<1>(result);

                reqVarInfo_t tmps = {EVALUATE, values, is_valid};
                requestedVars_.insert(
                                      std::pair<std::string, reqVarInfo_t>(namedIt->first, tmps));

            }

            // .. and constraints
            for(namedIt=constraints_.begin(); namedIt!=constraints_.end(); ++namedIt) {

                Expressions::Expr_t *constraint = namedIt->second;

                // find out which variables we need in order to evaluate the constraint
                variableDictionary_t variable_dictionary;
                getVariableDictionary(variable_dictionary,fn,constraint);

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
                                      std::pair<std::string, reqVarInfo_t>(namedIt->first, tmps));

            }
        } catch(SDDSParserException &e) {
            std::cout << "Evaluation of objective or constraint " << namedIt->first << " threw an exception ('" << e.what() << "' in " << e.where() << ")!" << std::endl;
            invalidBunch();
        } catch(OptPilotException &e) {
            std::cout << "Evaluation of objective or constraint " << namedIt->first << " threw an exception ('" << e.what() << "' in " << e.where() << ")!" << std::endl;
            invalidBunch();
        } catch(std::exception &e) {
            std::cout << "Evaluation of objective or constraint " << namedIt->first << " threw an exception ('" << e.what() << "')!" << std::endl;
            invalidBunch();
        } catch(...) {
            std::cout << "Evaluation of objective or constraint " << namedIt->first << " threw an exception!" << std::endl;
            invalidBunch();
        }

    }

    err = chdir(pwd_.c_str());
    if (err != 0) {
        std::cout << "Cannot chdir to "
                  << simulationDirName_.c_str() << std::endl;
    }
}

void OpalSimulation::getVariableDictionary(variableDictionary_t& dictionary,
                                           const std::string& filename,
                                           const Expressions::Expr_t* const expression) {

    std::set<std::string> req_vars = expression->getReqVars();

    // first check if required variables are design variables
    for (auto req_it = req_vars.begin(); req_it!=req_vars.end();) {
        auto it = userVariables_.find(*req_it);
        if (it==userVariables_.end()) { // not a design var
            ++req_it;
            continue;
        }
        double value = std::stod((*it).second);
        dictionary.insert(std::pair<std::string, double>(*req_it, value));
        req_it = req_vars.erase(req_it); // remove and update iterator to next
    }

    if(req_vars.empty()) return;

    // get remaining required variable values from the stat file
    boost::scoped_ptr<SDDSReader> sddsr(new SDDSReader(filename));
    sddsr->parseFile();

    for(std::string req_var : req_vars) {
        if(dictionary.count(req_var) != 0) continue;

        double value = 0.0;
        sddsr->getValue(-1 /*atTime*/, req_var, value);
        dictionary.insert(std::pair<std::string, double>(req_var, value));
    }
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
    namespace fs = boost::filesystem;
    try {
        int my_rank = 0;
        MPI_Comm_rank(comm_, &my_rank);
        if (my_rank == 0) {
            fs::path p(simulationDirName_.c_str());
            fs::remove_all(p);
        }
    } catch(fs::filesystem_error &ex) {
        std::cout << "Can't remove directory '" << simulationDirName_ << "', (" << ex.what() << ")" << std::endl;
    } catch(...) {
        std::cout << "Can't remove directory '" << simulationDirName_ << "'" << std::endl;
    }
}

void OpalSimulation::cleanUp(const std::vector<std::string>& keep) {
    namespace fs = boost::filesystem;

    if ( keep.empty() ) {
        // if empty we keep all files
        return;
    }

    try {
        int my_rank = 0;
        MPI_Comm_rank(comm_, &my_rank);
        if (my_rank == 0) {
            fs::path p(simulationDirName_.c_str());
            fs::directory_iterator it{p};
            while ( it != fs::directory_iterator{} ) {
                std::string extension = Util::toUpper(fs::extension(it->path().filename()));

                // remove .
                extension.erase(0, 1);

                auto result = std::find(keep.begin(), keep.end(), extension);

                if ( result == keep.end() ) {
                    fs::remove(it->path());
                }
                ++it;
            }
        }
    } catch(fs::filesystem_error &ex) {
        std::cout << "Can't remove file in directory '" << simulationDirName_
                  << "', (" << ex.what() << ")" << std::endl;
    } catch(...) {
        std::cout << "Can't remove file in directory '" << simulationDirName_
                  << "'" << std::endl;
    }
}