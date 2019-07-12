#include "Sample/SampleCmd.h"
#include "Sample/Sampler.h"
#include "Sample/OpalSample.h"
#include "Sample/RNGStream.h"

#include "Optimize/DVar.h"
#include "Optimize/Objective.h"
#include "Optimize/OpalSimulation.h"

#include "Attributes/Attributes.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include "Utility/IpplInfo.h"
#include "Utility/IpplTimings.h"
#include "Track/Track.h"

#include "Sample/SamplePilot.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"

#include "Comm/CommSplitter.h"
#include "Comm/Topology/NoCommTopology.h"
#include "Comm/Splitter/ManyMasterSplit.h"
#include "Comm/MasterGraph/SocialNetworkGraph.h"

#include "Expression/Parser/function.hpp"
#include "Expression/FromFile.h"
#include "Expression/SumErrSq.h"
#include "Expression/SDDSVariable.h"
#include "Expression/RadialPeak.h"
#include "Expression/MaxNormRadialPeak.h"
#include "Expression/NumberOfPeaks.h"
#include "Expression/SumErrSqRadialPeak.h"
#include "Expression/ProbeVariable.h"

#include <boost/filesystem.hpp>

#include <map>
#include <set>
#include <string>
#include <vector>

extern Inform *gmsg;

namespace {
    enum {
        INPUT,
        OUTPUT,
        OUTDIR,
        OBJECTIVES,
        STOREOBJECTIVES,
        DVARS,
        SAMPLINGS,
        NUMMASTERS,
        NUMCOWORKERS,
        TEMPLATEDIR,
        FIELDMAPDIR,
        DISTDIR,
        RASTER,
        SEED,
        KEEP,
        RESTART_FILE,
        RESTART_STEP,
        JSON_DUMP_FREQ,
        SIZE
    };
}

SampleCmd::SampleCmd():
    Action(SIZE, "SAMPLE",
           "The \"SAMPLE\" command initiates sampling.") {
    itsAttr[INPUT] = Attributes::makeString
        ("INPUT", "Path to input file");
    itsAttr[OUTPUT] = Attributes::makeString
        ("OUTPUT", "Name used in output file sample");
    itsAttr[OUTDIR] = Attributes::makeString
        ("OUTDIR", "Name of directory used to run and store sample output files");
    itsAttr[OBJECTIVES] = Attributes::makeStringArray
        ("OBJECTIVES", "List of expressions to evaluate and store");
    itsAttr[STOREOBJECTIVES] = Attributes::makeStringArray
        ("STOREOBJECTIVES", "List of stat variables to store");
    itsAttr[DVARS] = Attributes::makeStringArray
        ("DVARS", "List of sampling variables to be used");
    itsAttr[SAMPLINGS] = Attributes::makeStringArray
        ("SAMPLINGS", "List of sampling methods to be used");
    itsAttr[NUMMASTERS] = Attributes::makeReal
        ("NUM_MASTERS", "Number of master nodes");
    itsAttr[NUMCOWORKERS] = Attributes::makeReal
        ("NUM_COWORKERS", "Number processors per worker");
    itsAttr[TEMPLATEDIR] = Attributes::makeString
        ("TEMPLATEDIR", "Directory where templates are stored");
    itsAttr[FIELDMAPDIR] = Attributes::makeString
        ("FIELDMAPDIR", "Directory where field maps are stored");
    itsAttr[DISTDIR] = Attributes::makeString
        ("DISTDIR", "Directory where distributions are stored");
    itsAttr[RASTER] = Attributes::makeBool
        ("RASTER", "Scan full space given by design variables (default: true)", true);
    itsAttr[SEED] = Attributes::makeReal
        ("SEED", "Seed for global random number generator (default: 42)", 42);
    itsAttr[KEEP] = Attributes::makeStringArray
        ("KEEP", "List of files to keep for each simulation. (default: all files kept)");
    itsAttr[RESTART_FILE] = Attributes::makeString
        ("RESTART_FILE", "H5 file to restart the OPAL simulations from (optional)", "");
    itsAttr[RESTART_STEP] = Attributes::makeReal
        ("RESTART_STEP", "Restart from given H5 step (optional)",
         std::numeric_limits<int>::min());
    itsAttr[JSON_DUMP_FREQ] = Attributes::makeReal
        ("JSON_DUMP_FREQ", "Defines how often new individuals are appended to the final JSON file, "
         "i.e. every time JSON_DUMP_FREQ samples finished they are written (optional)",
         std::numeric_limits<size_t>::max());
    registerOwnership(AttributeHandler::COMMAND);
}

SampleCmd::SampleCmd(const std::string &name, SampleCmd *parent):
    Action(name, parent)
{ }

SampleCmd::~SampleCmd()
{ }

SampleCmd *SampleCmd::clone(const std::string &name) {
    return new SampleCmd(name, this);
}

void SampleCmd::execute() {

    namespace fs = boost::filesystem;

    auto opal = OpalData::getInstance();
    opal->setOptimizerFlag();

    fs::path inputfile(Attributes::getString(itsAttr[INPUT]));

    unsigned int seed = Attributes::getReal(itsAttr[SEED]);
    RNGStream::setGlobalSeed(seed);

    std::vector<std::string> objectivesstr  = Attributes::getStringArray(itsAttr[OBJECTIVES]);
    // FIXME Open issue #250 (https://gitlab.psi.ch/OPAL/src/issues/250)
    std::vector<std::string> storeobjstr  = Attributes::getStringArray(itsAttr[STOREOBJECTIVES]);
    std::vector<std::string> dvarsstr = Attributes::getStringArray(itsAttr[DVARS]);
    Expressions::Named_t objectives;
    DVarContainer_t dvars;

    std::vector<std::string> filesToKeep = Attributes::getStringArray(itsAttr[KEEP]);
    std::vector<std::string> sampling = Attributes::getStringArray(itsAttr[SAMPLINGS]);

    if ( sampling.size() != dvarsstr.size() )
        throw OpalException("SampleCmd::execute",
                            "Number of sampling methods != number of design variables.");

    typedef std::map< std::string, std::shared_ptr<SamplingMethod> > sampleMethods_t;
    sampleMethods_t sampleMethods;

    std::map<std::string, std::pair<double, double> > vars;

    for (std::string &name : dvarsstr) {
        name = Util::toUpper(name);
        Object *obj = opal->find(name);
        DVar* dvar = dynamic_cast<DVar*>(obj);
        if (dvar == nullptr) {
            throw OpalException("SampleCmd::execute",
                                "The sampling variable " + name + " is not known");

        }
        std::string var = dvar->getVariable();
        double lowerbound = dvar->getLowerBound();
        double upperbound = dvar->getUpperBound();

        auto ret = vars.insert(std::make_pair(var, std::make_pair(lowerbound, upperbound)));
        if (ret.second == false) {
            throw OpalException("SampleCmd::execute",
                                "There is already a design variable with the variable " + var + " defined");
        }

        DVar_t tmp = boost::make_tuple(var, lowerbound, upperbound);
        dvars.insert(namedDVar_t(name, tmp));
    }

    //////////////////////////////////////////////////////////////////////////

    functionDictionary_t funcs;
    client::function::type ff;
    ff = FromFile();
    funcs.insert(std::pair<std::string, client::function::type>
                 ("fromFile", ff));

    ff = SumErrSq();
    funcs.insert(std::pair<std::string, client::function::type>
                 ("sumErrSq", ff));

    ff = SDDSVariable();
    funcs.insert(std::pair<std::string, client::function::type>
                 ("sddsVariableAt", ff));

    ff = RadialPeak();
    funcs.insert(std::pair<std::string, client::function::type>
                 ("radialPeak", ff));

    ff = MaxNormRadialPeak();
    funcs.insert(std::pair<std::string, client::function::type>
                 ("maxNormRadialPeak", ff));

    ff = NumberOfPeaks();
    funcs.insert(std::pair<std::string, client::function::type>
            ("numberOfPeaks", ff));

    ff = SumErrSqRadialPeak();
    funcs.insert(std::pair<std::string, client::function::type>
                 ("sumErrSqRadialPeak", ff));

    ff = ProbeVariable();
    funcs.insert(std::pair<std::string, client::function::type>
                 ("probVariableWithID", ff));

    std::string fname = inputfile.stem().native();
    ff = sameSDDSVariable(fname);
    funcs.insert(std::pair<std::string, client::function::type>
                 ("statVariableAt", ff));

    //////////////////////////////////////////////////////////////////////////

    std::set<std::string> objExpressions; // check if all unique objective expressions
    for (std::string name: objectivesstr) {
        name = Util::toUpper(name);
        Object *obj = opal->find(name);
        Objective* objective = dynamic_cast<Objective*>(obj);
        if (objective == nullptr) {
            throw OpalException("OptimizeCmd::execute",
                                "The objective " + name + " is not known");

        }
        std::string expr = objective->getExpression();
        objectives.insert(Expressions::SingleNamed_t(
                   name, new Expressions::Expr_t(expr, funcs)));
        auto ret = objExpressions.insert(expr);
        if (ret.second == false) {
            throw OpalException("OptimizeCmd::execute",
                                "There is already a objective with the expression " + expr + " defined");
        }
    }

    bool raster = Attributes::getBool(itsAttr[RASTER]);
    size_t modulo = 1;
    unsigned int nSample = std::numeric_limits<unsigned int>::max();

    std::set<std::string> names; // check if all unique variables
    for (size_t i = 0; i < sampling.size(); ++i) {
        // corresponding sampling method
        OpalSample *s = OpalSample::find(sampling[i]);
        if (s == nullptr) {
            throw OpalException("SampleCmd::execute",
                                "Sampling method not found.");
        }

        std::string name = s->getVariable();

        if ( vars.find(name) == vars.end() ) {
            throw OpalException("SampleCmd::execute",
                                "Variable '" + name + "' not a DVAR.");
        }

        auto ret = names.insert(name);
        if (ret.second == false) {
            throw OpalException("SampleCmd::execute",
                                "There is already a sampling method with the variable " + name + " defined");
        }

        s->initialize(name,
                      vars[name].first,
                      vars[name].second,
                      modulo,
                      raster);

        if ( raster )
            modulo *= s->getSize();

        nSample = std::min(nSample, s->getSize());

        sampleMethods[name] = s->sampleMethod_m;
    }

    // Setup/Configuration
    //////////////////////////////////////////////////////////////////////////
    typedef OpalSimulation Sim_t;

    typedef CommSplitter< ManyMasterSplit< NoCommTopology > > Comm_t;
    typedef SocialNetworkGraph< NoCommTopology > SolPropagationGraph_t;

    typedef SamplePilot<Sampler, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

    //////////////////////////////////////////////////////////////////////////

    std::vector<std::string> arguments(opal->getArguments());
    std::vector<char*> argv;
    std::map<unsigned int, std::string> argumentMapper({
            {INPUT, "inputfile"},
            {OUTPUT, "outfile"},
            {OUTDIR, "outdir"},
            {NUMMASTERS, "num-masters"},
            {NUMCOWORKERS, "num-coworkers"},
            {RESTART_FILE, "restartfile"},
            {RESTART_STEP, "restartstep"},
            {JSON_DUMP_FREQ, "jsonDumpFreq"}
        });

    auto it = argumentMapper.end();
    for (unsigned int i = 0; i < SIZE; ++ i) {
        if ((it = argumentMapper.find(i)) != argumentMapper.end()) {
            std::string type = itsAttr[i].getType();
            if (type == "string") {
                if (Attributes::getString(itsAttr[i]) != "") {
                    std::string argument = "--" + (*it).second + "=" + Attributes::getString(itsAttr[i]);
                    arguments.push_back(argument);
                }
            } else if (type == "real") {
                if (itsAttr[i]) {
                    std::string val = std::to_string (Attributes::getReal(itsAttr[i]));
                    size_t last = val.find_last_not_of('0');
                    if (val[last] != '.') ++ last;
                    val.erase (last, std::string::npos );
                    std::string argument = "--" + (*it).second + "=" + val;
                    arguments.push_back(argument);
                }
            } else if (type == "logical") {
                if (itsAttr[i]) {
                    std::string argument = "--" + (*it).second + "=" + std::to_string(Attributes::getBool(itsAttr[i]));
                    arguments.push_back(argument);
                }
            }
        }
    }

    if ( raster )
        nSample = modulo;

    arguments.push_back("--nsamples=" + std::to_string(nSample));

    if (Attributes::getString(itsAttr[INPUT]) == "") {
        throw OpalException("SampleCmd::execute",
                            "The argument INPUT has to be provided");
    }

    if (Attributes::getString(itsAttr[OUTDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[OUTDIR]));
        if (dir.is_relative()) {
            fs::path path = fs::path(std::string(getenv("PWD")));
            path /= dir;
            dir = path;
        }

        if (!fs::exists(dir)) {
            fs::create_directory(dir);
        }
        std::string argument = "--simtmpdir=" + dir.native();
        arguments.push_back(argument);
    }

    if (Attributes::getString(itsAttr[TEMPLATEDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[TEMPLATEDIR]));
        if (dir.is_relative()) {
            fs::path path = fs::path(std::string(getenv("PWD")));
            path /= dir;
            dir = path;
        }

        std::string argument = "--templates=" + dir.native();
        arguments.push_back(argument);
    }

    if (Attributes::getString(itsAttr[FIELDMAPDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[FIELDMAPDIR]));
        if (dir.is_relative()) {
            fs::path path = fs::path(std::string(getenv("PWD")));
            path /= dir;
            dir = path;
        }

        setenv("FIELDMAPS", dir.c_str(), 1);
    }

    if (Attributes::getString(itsAttr[DISTDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[DISTDIR]));
        if (dir.is_relative()) {
            fs::path path = fs::path(std::string(getenv("PWD")));
            path /= dir;
            dir = path;
        }

        setenv("DISTRIBUTIONS", dir.c_str(), 1);
    }

    {
        std::string tmplFile = Attributes::getString(itsAttr[INPUT]);
        size_t pos = tmplFile.find_last_of("/");
        if(pos != std::string::npos)
            tmplFile = tmplFile.substr(pos+1);
        pos = tmplFile.find(".");
        tmplFile = tmplFile.substr(0,pos);
        tmplFile = Attributes::getString(itsAttr[TEMPLATEDIR]) + "/" + tmplFile + ".tmpl";

        std::ifstream infile(tmplFile.c_str());

        std::map<std::string, short> dvarCheck;
        auto itr = dvars.begin();
        for (; itr != dvars.end(); ++ itr) {
            dvarCheck.insert(std::make_pair(boost::get<0>(itr->second), 0));
        }

        while(infile.good()) {
            std::string line;
            std::getline(infile, line, '\n');

            //XXX doing the inverse would be better
            for(auto &check: dvarCheck) {
                size_t pos = line.find("_" + check.first + "_");
                if (pos != std::string::npos &&
                    dvarCheck.find(check.first) != dvarCheck.end()) {
                    dvarCheck.at(check.first) = 1;
                }
            }
        }
        infile.close();

        for (auto itr = dvarCheck.begin(); itr != dvarCheck.end(); ++ itr) {
            if (itr->second == 0) {
                throw OpalException("SampleCmd::execute()",
                                    "Couldn't find the design variable '" + itr->first + "' in '" + tmplFile + "'!");
            }
        }
    }

    *gmsg << endl;
    for (size_t i = 0; i < arguments.size(); ++ i) {
        argv.push_back(const_cast<char*>(arguments[i].c_str()));
        *gmsg << arguments[i] << " ";
    }
    *gmsg << endl;

    Inform *origGmsg = gmsg;
    gmsg = 0;
    try {
        CmdArguments_t args(new CmdArguments(argv.size(), &argv[0]));

        boost::shared_ptr<Comm_t>  comm(new Comm_t(args, MPI_COMM_WORLD));
        if (comm->isWorker())
            stashEnvironment();
        
        if ( comm->isOptimizer() ) {
            for (sampleMethods_t::iterator it = sampleMethods.begin();
                 it != sampleMethods.end(); ++it)
            {
                it->second->allocate(args, comm->getBundle());
            }
        }

        boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm, funcs, dvars,
                                                  objectives, sampleMethods,
                                                  storeobjstr, filesToKeep));
        if (comm->isWorker())
            popEnvironment();

    } catch (OptPilotException &e) {
        std::cout << "Exception caught: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -100);
    }
    gmsg = origGmsg;
}

void SampleCmd::stashEnvironment() {
    Ippl::stash();
    IpplTimings::stash();
    Track::stash();
    OpalData::stashInstance();
}

void SampleCmd::popEnvironment() {
    Ippl::pop();
    IpplTimings::pop();
    OpalData::popInstance();
    Track::pop();
}