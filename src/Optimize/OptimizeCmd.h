#ifndef OPAL_OptimizeCmd_HH
#define OPAL_OptimizeCmd_HH

#include "AbstractObjects/Action.h"

#include "Util/CmdArguments.h"
#include "Optimize/DVar.h"
#include "Expression/Expression.h"

#include <string>

// Class OptimizeCmd
// ------------------------------------------------------------------------
/// The OPTIMIZE command.

class OptimizeCmd: public Action {

public:

    /// Exemplar constructor.
    OptimizeCmd();

    virtual ~OptimizeCmd();

    /// Make clone.
    virtual OptimizeCmd *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    OptimizeCmd(const OptimizeCmd &)    = delete;
    void operator=(const OptimizeCmd &) = delete;

    // Clone constructor.
    OptimizeCmd(const std::string &name, OptimizeCmd *parent);

    void stashEnvironment();
    void popEnvironment();
    
    enum CrossOver {
        Blend = 0,
        NaiveOnePoint,
        NaiveUniform,
        SimulatedBinary
    };
    
    CrossOver crossoverSelection(std::string crossover);
    
    enum Mutation {
        IndependentBit = 10,
        OneBit = 20
    };
    
    Mutation mutationSelection(std::string mutation);
    
    void run(const CmdArguments_t& args,
             const functionDictionary_t& funcs,
             const DVarContainer_t& dvars,
             const Expressions::Named_t& objectives,
             const Expressions::Named_t& constraints,
             const std::vector<std::string>&);
};

#endif