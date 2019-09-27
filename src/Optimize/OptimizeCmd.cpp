#include "Optimize/OptimizeCmd.h"
#include "Optimize/DVar.h"
#include "Optimize/Objective.h"
#include "Optimize/Constraint.h"
#include "Optimize/OpalSimulation.h"

#include "Attributes/Attributes.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"

//#include "Utility/Inform.h"
#include "Utility/IpplInfo.h"
#include "Utility/IpplTimings.h"
#include "Track/Track.h"

#include "Pilot/Pilot.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"
#include "Util/OpalInputFileParser.h"

#include "Optimizer/EA/FixedPisaNsga2.h"
#include "Optimizer/EA/BlendCrossover.h"
#include "Optimizer/EA/IndependentBitMutation.h"

#include "Comm/CommSplitter.h"
#include "Comm/Topology/NoCommTopology.h"
#include "Comm/Splitter/ManyMasterSplit.h"
#include "Comm/MasterGraph/SocialNetworkGraph.h"

#include "Expression/Parser/function.hpp"
#include "Expression/FromFile.h"
#include "Expression/SumErrSq.h"
#include "Expression/SDDSVariable.h"
#include "Expression/RadialPeak.h"
#include "Expression/SumErrSqRadialPeak.h"
#include "Expression/ProbeVariable.h"

#include <boost/filesystem.hpp>

extern Inform *gmsg;

namespace {
    enum {
        INPUT,
        OUTPUT,
        OUTDIR,
        OBJECTIVES,
        DVARS,
        CONSTRAINTS,
        INITIALPOPULATION,
        NUMMASTERS,
        NUMCOWORKERS,
        SELECTOR,
        DUMPDAT,
        DUMPFREQ,
        NUMINDGEN,
        MAXGENERATIONS,
        EPSILON,
        EXPECTEDHYPERVOL,
        CONVHVOLPROG,
        ONEPILOTCONVERGE,
        SOLSYNCH,
        GENEMUTATIONPROBABILITY,
        MUTATIONPROBABILITY,
        RECOMBINATIONPROBABILITY,
        SIMBINCROSSOVERNU,
        SIMTMPDIR,
        TEMPLATEDIR,
        FIELDMAPDIR,
        KEEP,
        SIZE
    };
}

OptimizeCmd::OptimizeCmd():
    Action(SIZE, "OPTIMIZE",
           "The \"OPTIMIZE\" command initiates optimization.") {
    itsAttr[INPUT] = Attributes::makeString
        ("INPUT", "Path to input file");
    itsAttr[OUTPUT] = Attributes::makeString
        ("OUTPUT", "Name used in output file generation");
    itsAttr[OUTDIR] = Attributes::makeString
        ("OUTDIR", "Name of directory used to store generation output files");
    itsAttr[OBJECTIVES] = Attributes::makeStringArray
        ("OBJECTIVES", "List of objectives to be used");
    itsAttr[DVARS] = Attributes::makeStringArray
        ("DVARS", "List of optimization variables to be used");
    itsAttr[CONSTRAINTS] = Attributes::makeStringArray
        ("CONSTRAINTS", "List of constraints to be used");
    itsAttr[INITIALPOPULATION] = Attributes::makeReal
        ("INITIALPOPULATION", "Size of the initial population");
    itsAttr[NUMMASTERS] = Attributes::makeReal
        ("NUM_MASTERS", "Number of master nodes");
    itsAttr[NUMCOWORKERS] = Attributes::makeReal
        ("NUM_COWORKERS", "Number processors per worker");
    itsAttr[SELECTOR] = Attributes::makeString
        ("SELECTOR", "Path of the selector (PISA only)");
    itsAttr[DUMPDAT] = Attributes::makeReal
        ("DUMP_DAT", "Dump old generation data format with frequency (PISA only)");
    itsAttr[DUMPFREQ] = Attributes::makeReal
        ("DUMP_FREQ", "Dump old generation data format with frequency (PISA only)");
    itsAttr[NUMINDGEN] = Attributes::makeReal
        ("NUM_IND_GEN", "Number of individuals in a generation (PISA only)");
    itsAttr[MAXGENERATIONS] = Attributes::makeReal
        ("MAXGENERATIONS", "Number of generations to run");
    itsAttr[EPSILON] = Attributes::makeReal
        ("EPSILON", "Tolerance of hypervolume criteria");
    itsAttr[EXPECTEDHYPERVOL] = Attributes::makeReal
        ("EXPECTED_HYPERVOL", "The reference hypervolume");
    itsAttr[CONVHVOLPROG] = Attributes::makeReal
        ("CONV_HVOL_PROG", "converge if change in hypervolume is smaller");
    itsAttr[ONEPILOTCONVERGE] = Attributes::makeBool
        ("ONE_PILOT_CONVERGE", "");
    itsAttr[SOLSYNCH] = Attributes::makeReal
        ("SOL_SYNCH", "Solution exchange frequency");
    itsAttr[GENEMUTATIONPROBABILITY] = Attributes::makeReal
        ("GENE_MUTATION_PROBABILITY", "Mutation probability of individual gene, default: 0.5");
    itsAttr[MUTATIONPROBABILITY] = Attributes::makeReal
        ("MUTATION_PROBABILITY", "Mutation probability of genom, default: 0.5");
    itsAttr[RECOMBINATIONPROBABILITY] = Attributes::makeReal
        ("RECOMBINATION_PROBABILITY", "Probability for genes to recombine, default: 0.5");
    itsAttr[SIMBINCROSSOVERNU] = Attributes::makeReal
        ("SIMBIN_CROSSOVER_NU", "Simulated binary crossover");
    itsAttr[SIMTMPDIR] = Attributes::makeString
        ("SIMTMPDIR", "Directory where simulations are run");
    itsAttr[TEMPLATEDIR] = Attributes::makeString
        ("TEMPLATEDIR", "Directory where templates are stored");
    itsAttr[FIELDMAPDIR] = Attributes::makeString
        ("FIELDMAPDIR", "Directory where field maps are stored");
    itsAttr[KEEP] = Attributes::makeStringArray
        ("KEEP", "List of files to keep for each simulation. (default: all files kept)");


    registerOwnership(AttributeHandler::COMMAND);
}

OptimizeCmd::OptimizeCmd(const std::string &name, OptimizeCmd *parent):
    Action(name, parent)
{ }

OptimizeCmd::~OptimizeCmd()
{ }

OptimizeCmd *OptimizeCmd::clone(const std::string &name) {
    return new OptimizeCmd(name, this);
}

void OptimizeCmd::execute() {
    namespace fs = boost::filesystem;

    auto opal = OpalData::getInstance();
    fs::path inputfile(Attributes::getString(itsAttr[INPUT]));

    std::vector<std::string> filesToKeep    = Attributes::getStringArray(itsAttr[KEEP]);
    std::vector<std::string> dvarsstr = Attributes::getStringArray(itsAttr[DVARS]);
    std::vector<std::string> objectivesstr = Attributes::getStringArray(itsAttr[OBJECTIVES]);
    std::vector<std::string> constraintsstr = Attributes::getStringArray(itsAttr[CONSTRAINTS]);
    DVarContainer_t dvars;
    Expressions::Named_t objectives;
    Expressions::Named_t constraints;

    // Setup/Configuration
    //////////////////////////////////////////////////////////////////////////
    typedef OpalInputFileParser Input_t;
    typedef OpalSimulation Sim_t;

    typedef FixedPisaNsga2< BlendCrossover, IndependentBitMutation > Opt_t;

    typedef CommSplitter< ManyMasterSplit< NoCommTopology > > Comm_t;
    typedef SocialNetworkGraph< NoCommTopology > SolPropagationGraph_t;

    typedef Pilot<Input_t, Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

    // prepare function dictionary and add all available functions in
    // expressions
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

    std::vector<std::string> arguments(opal->getArguments());
    std::vector<char*> argv;
    std::map<unsigned int, std::string> argumentMapper({
            {INPUT, "inputfile"},
            {OUTPUT, "outfile"},
            {OUTDIR, "outdir"},
            {INITIALPOPULATION, "initialPopulation"},
            {NUMMASTERS, "num-masters"},
            {NUMCOWORKERS, "num-coworkers"},
            {SELECTOR, "selector"},
            {DUMPDAT, "dump-dat"},
            {DUMPFREQ, "dump-freq"},
            {NUMINDGEN, "num-ind-gen"},
            {MAXGENERATIONS, "maxGenerations"},
            {EPSILON, "epsilon"},
            {EXPECTEDHYPERVOL, "expected-hypervol"},
            {CONVHVOLPROG, "conv-hvol-prog"},
            {ONEPILOTCONVERGE, "one-pilot-converge"},
            {SOLSYNCH, "sol-synch"},
            {GENEMUTATIONPROBABILITY, "gene-mutation-probability"},
            {MUTATIONPROBABILITY, "mutation-probability"},
            {RECOMBINATIONPROBABILITY, "recombination-probability"},
            {SIMBINCROSSOVERNU, "simbin-crossover-nu"}
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
    if (Attributes::getString(itsAttr[INPUT]) == "") {
        throw OpalException("OptimizeCmd::execute",
                            "The argument INPUT has to be provided");
    }
    if (Attributes::getReal(itsAttr[INITIALPOPULATION]) <= 0) {
        throw OpalException("OptimizeCmd::execute",
                            "The argument INITIALPOPULATION has to be provided");
    }
    if (Attributes::getReal(itsAttr[MAXGENERATIONS]) <= 0) {
        throw OpalException("OptimizeCmd::execute",
                            "The argument MAXGENERATIONS has to be provided");
    }

    if (Attributes::getString(itsAttr[SIMTMPDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[SIMTMPDIR]));
        if (dir.is_relative()) {
            fs::path path = fs::path(std::string(getenv("PWD")));
            path /= dir;
            dir = path;
        }

        *gmsg << dir.native() << endl;
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

    *gmsg << endl;
    for (size_t i = 0; i < arguments.size(); ++ i) {
        argv.push_back(const_cast<char*>(arguments[i].c_str()));
        *gmsg << arguments[i] << " ";
    }
    *gmsg << endl;

    for (const std::string &name: dvarsstr) {
        Object *obj = opal->find(name);
        DVar* dvar = static_cast<DVar*>(obj);
        std::string var = dvar->getVariable();
        double lowerbound = dvar->getLowerBound();
        double upperbound = dvar->getUpperBound();

        DVar_t tmp = boost::make_tuple(var, lowerbound, upperbound);
        dvars.insert(namedDVar_t(name, tmp));
    }
    for (const std::string &name: objectivesstr) {
        Object *obj = opal->find(name);
        Objective* objective = static_cast<Objective*>(obj);
        std::string expr = objective->getExpression();
        objectives.insert(Expressions::SingleNamed_t(
                   name, new Expressions::Expr_t(expr, funcs)));
    }
    for (const std::string &name: constraintsstr) {
        Object *obj = opal->find(name);
        Constraint* constraint = static_cast<Constraint*>(obj);
        std::string expr = constraint->getExpression();
        constraints.insert(Expressions::SingleNamed_t(
                    name, new Expressions::Expr_t(expr, funcs)));
    }

    Inform *origGmsg = gmsg;
    gmsg = 0;
    stashEnvironment();
    try {
        CmdArguments_t args(new CmdArguments(argv.size(), &argv[0]));

        boost::shared_ptr<Comm_t>  comm(new Comm_t(args, MPI_COMM_WORLD));
        boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm, funcs, dvars, objectives, constraints, filesToKeep));

    } catch (OptPilotException &e) {
        std::cout << "Exception caught: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -100);
    }
    popEnvironment();
    gmsg = origGmsg;
}

void OptimizeCmd::stashEnvironment() {
    Ippl::stash();
    IpplTimings::stash();
    Track::stash();
    OpalData::stashInstance();
}

void OptimizeCmd::popEnvironment() {
    Ippl::pop();
    IpplTimings::pop();
    OpalData::popInstance();
    Track::pop();
}