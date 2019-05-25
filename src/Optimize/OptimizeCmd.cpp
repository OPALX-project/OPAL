#include "Optimize/OptimizeCmd.h"
#include "Optimize/Objective.h"
#include "Optimize/Constraint.h"
#include "Optimize/OpalSimulation.h"

#include "Attributes/Attributes.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

//#include "Utility/Inform.h"
#include "Utility/IpplInfo.h"
#include "Utility/IpplTimings.h"
#include "Track/Track.h"

#include "Pilot/Pilot.h"
#include "Util/OptPilotException.h"

#include "Optimizer/EA/FixedPisaNsga2.h"
#include "Optimizer/EA/BlendCrossover.h"
#include "Optimizer/EA/SimulatedBinaryCrossover.h"
#include "Optimizer/EA/NaiveOnePointCrossover.h"
#include "Optimizer/EA/NaiveUniformCrossover.h"
#include "Optimizer/EA/IndependentBitMutation.h"
#include "Optimizer/EA/OneBitMutation.h"

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
        DVARS,
        CONSTRAINTS,
        INITIALPOPULATION,
        STARTPOPULATION,
        NUMMASTERS,
        NUMCOWORKERS,
        DUMPDAT,
        DUMPFREQ,
        DUMPOFFSPRING,
        NUMINDGEN,
        MAXGENERATIONS,
        EPSILON,
        EXPECTEDHYPERVOL,
        HYPERVOLREFERENCE,
        CONVHVOLPROG,
        ONEPILOTCONVERGE,
        SOLSYNCH,
        GENEMUTATIONPROBABILITY,
        MUTATIONPROBABILITY,
        RECOMBINATIONPROBABILITY,
        SIMBINCROSSOVERNU,
        INITIALOPTIMIZATION,
        BIRTHCONTROL,
        SIMTMPDIR,
        TEMPLATEDIR,
        FIELDMAPDIR,
        DISTDIR,
        CROSSOVER,
        MUTATION,
        RESTART_FILE,
        RESTART_STEP,
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
    itsAttr[STARTPOPULATION] = Attributes::makeString
        ("STARTPOPULATION", "Generation file (JSON format) to be started from (optional)", "");
    itsAttr[NUMMASTERS] = Attributes::makeReal
        ("NUM_MASTERS", "Number of master nodes");
    itsAttr[NUMCOWORKERS] = Attributes::makeReal
        ("NUM_COWORKERS", "Number processors per worker");
    itsAttr[DUMPDAT] = Attributes::makeReal
        ("DUMP_DAT", "Dump old generation data format with frequency (PISA only)");
    itsAttr[DUMPFREQ] = Attributes::makeReal
        ("DUMP_FREQ", "Dump generation data with frequency (PISA only)");
    itsAttr[DUMPOFFSPRING] = Attributes::makeBool
        ("DUMP_OFFSPRING", "Dump offspring (instead of parent population), default: true");
    itsAttr[NUMINDGEN] = Attributes::makeReal
        ("NUM_IND_GEN", "Number of individuals in a generation (PISA only)");
    itsAttr[MAXGENERATIONS] = Attributes::makeReal
        ("MAXGENERATIONS", "Number of generations to run");
    itsAttr[EPSILON] = Attributes::makeReal
        ("EPSILON", "Tolerance of hypervolume criteria, default 0.001");
    itsAttr[EXPECTEDHYPERVOL] = Attributes::makeReal
        ("EXPECTED_HYPERVOL", "The reference hypervolume, default 0");
    itsAttr[HYPERVOLREFERENCE] = Attributes::makeRealArray
        ("HYPERVOLREFERENCE", "The reference point (real array) for the hypervolume, default empty (origin)");
    itsAttr[CONVHVOLPROG] = Attributes::makeReal
        ("CONV_HVOL_PROG", "Converge if change in hypervolume is smaller, default 0");
    itsAttr[ONEPILOTCONVERGE] = Attributes::makeBool
        ("ONE_PILOT_CONVERGE", "default false");
    itsAttr[SOLSYNCH] = Attributes::makeReal
        ("SOL_SYNCH", "Solution exchange frequency, default 0");
    itsAttr[GENEMUTATIONPROBABILITY] = Attributes::makeReal
        ("GENE_MUTATION_PROBABILITY", "Mutation probability of individual gene, default: 0.5");
    itsAttr[MUTATIONPROBABILITY] = Attributes::makeReal
        ("MUTATION_PROBABILITY", "Mutation probability of genome, default: 0.5");
    itsAttr[RECOMBINATIONPROBABILITY] = Attributes::makeReal
        ("RECOMBINATION_PROBABILITY", "Probability for genes to recombine, default: 0.5");
    itsAttr[SIMBINCROSSOVERNU] = Attributes::makeReal
        ("SIMBIN_CROSSOVER_NU", "Simulated binary crossover, default: 2.0");
    itsAttr[INITIALOPTIMIZATION] = Attributes::makeBool
        ("INITIAL_OPTIMIZATION", "Optimize speed of initial generation, default: false");
    itsAttr[BIRTHCONTROL] = Attributes::makeBool
        ("BIRTH_CONTROL", "Enforce strict population sizes (or flexible to keep workers busy), default: false");
    itsAttr[SIMTMPDIR] = Attributes::makeString
        ("SIMTMPDIR", "Directory where simulations are run");
    itsAttr[TEMPLATEDIR] = Attributes::makeString
        ("TEMPLATEDIR", "Directory where templates are stored");
    itsAttr[FIELDMAPDIR] = Attributes::makeString
        ("FIELDMAPDIR", "Directory where field maps are stored");
    itsAttr[DISTDIR] = Attributes::makeString
        ("DISTDIR", "Directory where distributions are stored", "");
    itsAttr[CROSSOVER] = Attributes::makeString
        ("CROSSOVER", "Type of cross over (default: Blend)", "Blend");
    itsAttr[MUTATION] = Attributes::makeString
        ("MUTATION", "Type of bit mutation (default: IndependentBit)", "IndependentBit");
    itsAttr[RESTART_FILE] = Attributes::makeString
        ("RESTART_FILE", "H5 file to restart the OPAL simulations from (optional)", "");
    itsAttr[RESTART_STEP] = Attributes::makeReal
        ("RESTART_STEP", "Restart from given H5 step (optional)",
         std::numeric_limits<int>::min());
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
    opal->setOptimizerFlag();

    fs::path inputfile(Attributes::getString(itsAttr[INPUT]));

    std::vector<std::string> dvarsstr       = Attributes::getStringArray(itsAttr[DVARS]);
    std::vector<std::string> objectivesstr  = Attributes::getStringArray(itsAttr[OBJECTIVES]);
    std::vector<std::string> constraintsstr = Attributes::getStringArray(itsAttr[CONSTRAINTS]);
    DVarContainer_t dvars;
    Expressions::Named_t objectives;
    Expressions::Named_t constraints;

    // Setup/Configuration
    //////////////////////////////////////////////////////////////////////////

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

    std::vector<std::string> arguments(opal->getArguments());
    std::vector<char*> argv;
    std::map<unsigned int, std::string> argumentMapper({
            {INPUT, "inputfile"},
            {OUTPUT, "outfile"},
            {OUTDIR, "outdir"},
            {INITIALPOPULATION, "initialPopulation"},
            {STARTPOPULATION, "start-population"},
            {NUMMASTERS, "num-masters"},
            {NUMCOWORKERS, "num-coworkers"},
            {DUMPDAT, "dump-dat"},
            {DUMPFREQ, "dump-freq"},
            {DUMPOFFSPRING, "dump-offspring"},
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
            {SIMBINCROSSOVERNU, "simbin-crossover-nu"},
            {INITIALOPTIMIZATION, "initial-optimization"},
            {BIRTHCONTROL, "birth-control"},
            {RESTART_FILE, "restartfile"},
            {RESTART_STEP, "restartstep"}
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
    // sanity checks
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
    if (Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE]).empty() == false &&
        Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE]).size()  != objectivesstr.size()) {
        throw OpalException("OptimizeCmd::execute",
                            "The hypervolume reference point should have the same dimension as the objectives");
    }
    if (Attributes::getString(itsAttr[STARTPOPULATION])     != "" &&
        Attributes::getBool(  itsAttr[INITIALOPTIMIZATION]) == true) {
        throw OpalException("OptimizeCmd::execute",
                            "No INITIAL_OPTIMIZATION possible when reading initial population from file (STARTPOPULATION)");
    }
    if (Attributes::getBool(itsAttr[BIRTHCONTROL])        == true &&
        Attributes::getBool(itsAttr[INITIALOPTIMIZATION]) == true) {
        throw OpalException("OptimizeCmd::execute",
                            "No INITIAL_OPTIMIZATION possible with BIRTH_CONTROL");
    }

    if (Attributes::getString(itsAttr[SIMTMPDIR]) != "") {
        fs::path dir(Attributes::getString(itsAttr[SIMTMPDIR]));
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


    *gmsg << endl;
    for (size_t i = 0; i < arguments.size(); ++ i) {
        argv.push_back(const_cast<char*>(arguments[i].c_str()));
        *gmsg << arguments[i] << " ";
    }
    *gmsg << endl;

    std::set<std::string> vars; // check if all unique vars
    for (const std::string &name: dvarsstr) {
        Object *obj = opal->find(name);
        DVar* dvar = dynamic_cast<DVar*>(obj);
        if (dvar == nullptr) {
            throw OpalException("OptimizeCmd::execute",
                                "The design variable " + name + " is not known");

        }
        std::string var = dvar->getVariable();
        double lowerbound = dvar->getLowerBound();
        double upperbound = dvar->getUpperBound();

        DVar_t tmp = boost::make_tuple(var, lowerbound, upperbound);
        dvars.insert(namedDVar_t(name, tmp));
        auto ret = vars.insert(var);
        if (ret.second == false) {
            throw OpalException("OptimizeCmd::execute",
                                "There is already a design variable with the variable " + var + " defined");
        }
    }
    std::set<std::string> objExpressions; // check if all unique objective expressions
    for (const std::string &name: objectivesstr) {
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
    std::set<std::string> constraintExpressions; // check if all unique constraint expressions
    for (const std::string &name: constraintsstr) {
        Object *obj = opal->find(name);
        Constraint* constraint = dynamic_cast<Constraint*>(obj);
        if (constraint == nullptr) {
            throw OpalException("OptimizeCmd::execute",
                                "The constraint " + name + " is not known");

        }
        std::string expr = constraint->getExpression();
        constraints.insert(Expressions::SingleNamed_t(
                    name, new Expressions::Expr_t(expr, funcs)));
        auto ret = constraintExpressions.insert(expr);
        if (ret.second == false) {
            throw OpalException("OptimizeCmd::execute",
                                "There is already a constraint with the expression " + expr + " defined");
        }
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
        for (auto itr = dvars.begin(); itr != dvars.end(); ++ itr) {
            dvarCheck.insert(std::make_pair(boost::get<0>(itr->second), 0));
        }

        while(infile.good()) {
            std::string line;
            std::getline(infile, line, '\n');

            //XXX doing the inverse would be better
            for(auto &check: dvarCheck) {
                pos = line.find("_" + check.first + "_");
                if (pos != std::string::npos &&
                    dvarCheck.find(check.first) != dvarCheck.end()) {
                    dvarCheck.at(check.first) = 1;
                }
            }
        }
        infile.close();

        for (auto itr = dvarCheck.begin(); itr != dvarCheck.end(); ++ itr) {
            if (itr->second == 0) {
                throw OpalException("OptimizeCmd::execute()",
                                    "Couldn't find the design variable '" + itr->first + "' in '" + tmplFile + "'!");
            }
        }
    }

    Inform *origGmsg = gmsg;
    gmsg = 0;
    try {
        CmdArguments_t args(new CmdArguments(argv.size(), &argv[0]));

        this->run(args, funcs, dvars, objectives, constraints);

    } catch (OptPilotException &e) {
        std::cout << "Exception caught: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -100);
    }
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

OptimizeCmd::CrossOver OptimizeCmd::crossoverSelection(std::string crossover) {
    crossover = Util::toUpper(crossover);

    std::map<std::string, CrossOver> map;
    map["BLEND"] = CrossOver::Blend;
    map["NAIVEONEPOINT"] = CrossOver::NaiveOnePoint;
    map["NAIVEUNIFORM"] = CrossOver::NaiveUniform;
    map["SIMULATEDBINARY"] = CrossOver::SimulatedBinary;

    CrossOver co = CrossOver::Blend;

    switch ( map[crossover] ) {
        case CrossOver::Blend:
            break;
        case CrossOver::NaiveOnePoint:
            co = CrossOver::NaiveOnePoint;
            break;
        case CrossOver::NaiveUniform:
            co = CrossOver::NaiveUniform;
            break;
        case CrossOver::SimulatedBinary:
            co = CrossOver::SimulatedBinary;
            break;
        default:
            throw OpalException("OptimizeCmd::crossoverSelection",
                                "No cross over '" + crossover + "' supported.");
    }

    return co;
}

OptimizeCmd::Mutation OptimizeCmd::mutationSelection(std::string mutation) {
    mutation = Util::toUpper(mutation);

    std::map<std::string, Mutation> map;
    map["INDEPENDENTBIT"] = Mutation::IndependentBit;
    map["ONEBIT"] = Mutation::OneBit;

    Mutation mut = Mutation::IndependentBit;

    switch ( map[mutation] ) {
        case Mutation::IndependentBit:
            break;
        case Mutation::OneBit:
            mut = Mutation::OneBit;
            break;
        default:
            throw OpalException("OptimizeCmd::mutationSelection",
                                "No mutation '" + mutation + "' supported.");
    }

    return mut;
}

void OptimizeCmd::run(const CmdArguments_t& args,
                      const functionDictionary_t& funcs,
                      const DVarContainer_t& dvars,
                      const Expressions::Named_t& objectives,
                      const Expressions::Named_t& constraints)
{
    typedef OpalSimulation Sim_t;

    typedef CommSplitter< ManyMasterSplit< NoCommTopology > > Comm_t;
    typedef SocialNetworkGraph< NoCommTopology > SolPropagationGraph_t;

    boost::shared_ptr<Comm_t>  comm(new Comm_t(args, MPI_COMM_WORLD));
    if (comm->isWorker())
        stashEnvironment();

    CrossOver crossover = this->crossoverSelection(Attributes::getString(itsAttr[CROSSOVER]));
    Mutation mutation = this->mutationSelection(Attributes::getString(itsAttr[MUTATION]));

    switch ( crossover + mutation ) {
        case CrossOver::Blend + Mutation::IndependentBit:
        {
            typedef FixedPisaNsga2< BlendCrossover, IndependentBitMutation > Opt_t;
            typedef Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

            boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm,
                                              funcs, dvars,
                                              objectives, constraints,
                                              Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE])));
            break;
        }
        case CrossOver::Blend + Mutation::OneBit:
        {
            typedef FixedPisaNsga2< BlendCrossover, OneBitMutation > Opt_t;
            typedef Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

            boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm,
                                              funcs, dvars,
                                              objectives, constraints,
                                              Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE])));
            break;
        }
        case CrossOver::NaiveOnePoint + Mutation::IndependentBit:
        {
            typedef FixedPisaNsga2< NaiveOnePointCrossover, IndependentBitMutation > Opt_t;
            typedef Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

            boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm,
                                              funcs, dvars,
                                              objectives, constraints,
                                              Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE])));
            break;
        }
        case CrossOver::NaiveOnePoint + Mutation::OneBit:
        {
            typedef FixedPisaNsga2< NaiveOnePointCrossover, OneBitMutation > Opt_t;
            typedef Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

            boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm,
                                              funcs, dvars,
                                              objectives, constraints,
                                              Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE])));
            break;
        }
        case CrossOver::NaiveUniform + Mutation::IndependentBit:
        {
            typedef FixedPisaNsga2< NaiveUniformCrossover, IndependentBitMutation > Opt_t;
            typedef Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

            boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm,
                                              funcs, dvars,
                                              objectives, constraints,
                                              Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE])));
            break;
        }
        case CrossOver::NaiveUniform + Mutation::OneBit:
        {
            typedef FixedPisaNsga2< NaiveUniformCrossover, OneBitMutation > Opt_t;
            typedef Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

            boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm,
                                              funcs, dvars,
                                              objectives, constraints,
                                              Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE])));
            break;
        }
        case CrossOver::SimulatedBinary + Mutation::IndependentBit:
        {
            typedef FixedPisaNsga2< SimulatedBinaryCrossover, IndependentBitMutation > Opt_t;
            typedef Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

            boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm,
                                              funcs, dvars,
                                              objectives, constraints,
                                              Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE])));
            break;
        }
        case CrossOver::SimulatedBinary + Mutation::OneBit:
        {
            typedef FixedPisaNsga2< SimulatedBinaryCrossover, OneBitMutation > Opt_t;
            typedef Pilot<Opt_t, Sim_t, SolPropagationGraph_t, Comm_t> pilot_t;

            boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm,
                                              funcs, dvars,
                                              objectives, constraints,
                                              Attributes::getRealArray(itsAttr[HYPERVOLREFERENCE])));
            break;
        }
        default:
            throw OpalException("OptimizeCmd::run",
                                "No such cross over and mutation combination supported.");
    }

    if (comm->isWorker())
        popEnvironment();
}