//
// Class OptimizeCmd
//   The OptimizeCmd definition.
//   A OptimizeCmd definition is used to parse the parametes for the optimizer.
//
// Copyright (c) 2017, Christof Metzger-Kraus
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPAL_OptimizeCmd_HH
#define OPAL_OptimizeCmd_HH

#include "AbstractObjects/Action.h"

#include "Util/CmdArguments.h"
#include "Optimize/DVar.h"
#include "Expression/Expression.h"

#include <string>

class OptimizeCmd: public Action {

public:

    /// Exemplar constructor.
    OptimizeCmd();

    virtual ~OptimizeCmd();

    /// Make clone.
    virtual OptimizeCmd *clone(const std::string& name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    OptimizeCmd(const OptimizeCmd&)    = delete;
    void operator=(const OptimizeCmd&) = delete;

    // Clone constructor.
    OptimizeCmd(const std::string &name, OptimizeCmd* parent);

    void stashEnvironment();
    void popEnvironment();

    enum class CrossOver {
        Blend,
        NaiveOnePoint,
        NaiveUniform,
        SimulatedBinary
    };
    CrossOver crossoverSelection(const std::string& crossover);

    enum class Mutation {
        IndependentBit,
        OneBit
    };
    Mutation mutationSelection(const std::string& mutation);

    void run(const CmdArguments_t& args,
             const functionDictionary_t& funcs,
             const DVarContainer_t& dvars,
             const Expressions::Named_t& objectives,
             const Expressions::Named_t& constraints);
};

#endif
