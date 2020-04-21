//
// Class Expression
//   Expression to be evaluated in the framework.
//
//   @see GlobalFunctions.h
//
//   This class uses the Boost Spirit parser to parse and evaluate string
//   expressions (objectives or constraints).
//   Custom functions called in the expression should be registered by the
//   driver. A collection of C math default functions is always included.
//   For constraints the operator type can be queried.
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
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
#ifndef __EXPRESSION_H__
#define __EXPRESSION_H__

#include <map>
#include <set>
#include <string>

#include "Util/Types.h"
#include "Util/OptPilotException.h"

#include "Expression/Parser/expression.hpp"
#include "Expression/Parser/evaluator.hpp"
#include "Expression/Parser/requirements.hpp"
#include "Expression/Parser/skipper.hpp"
#include "Expression/Parser/function.hpp"

#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/variant/get.hpp>
#include <boost/variant/variant.hpp>
#include "boost/smart_ptr.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/algorithm/string.hpp"


typedef std::map<std::string, double> variableDictionary_t;
typedef std::map<std::string, client::function::type> functionDictionary_t;


class Expression;
namespace Expressions {

    /// type of an expression
    typedef Expression Expr_t;

    // result of an evaluated expression
    typedef boost::tuple<double, bool> Result_t;
    enum Result_tIdx{
        VALUE,
        IS_VALID
    };

    //FIXME: this actually should be a map of type (name, fusion::vector<...>)
    /// type of an expressions with a name
    typedef std::map<std::string, Expressions::Expr_t*> Named_t;

    //XXX name of one single expression
    typedef std::pair<std::string, Expressions::Expr_t*> SingleNamed_t;

    /// distinguish different constraints
    enum OperatorType_t {
        NONE,
        EQ,             // ==
        NOT_EQ,         // !=
        INEQ_LHS,       // <
        INEQ_LHS_EQ,    // <=
        INEQ_RHS,       // >
        INEQ_RHS_EQ     // >=
    };

}

#include "Expression/GlobalFunctions.h"

class Expression {

public:

    Expression()
    {}

    Expression(std::string expr)
        : expr_(expr)
    {
        determineConstrOperator();
        functionDictionary_t global_funcs = GlobalFunctions::get();
        known_expr_funcs_.insert(global_funcs.begin(), global_funcs.end());
        parse();
    }

    Expression(std::string expr, functionDictionary_t known_expr_funcs)
        : expr_(expr)
        , known_expr_funcs_(known_expr_funcs)
    {
        determineConstrOperator();
        functionDictionary_t global_funcs = GlobalFunctions::get();
        known_expr_funcs_.insert(global_funcs.begin(), global_funcs.end());
        parse();
    }

    virtual ~Expression()
    {}


    std::set<std::string> getReqVars()  const { return vars_;  }
    std::set<std::string> getReqFuncs() const { return funcs_; }
    std::string toString()              const { return expr_;  }
    functionDictionary_t getRegFuncs()  const { return known_expr_funcs_; }

    /// get operator type present (if expression is constraint)
    Expressions::OperatorType_t getOpType() const { return type_; }

    /// evaluate an expression given a value dictionary of free variables
    Expressions::Result_t evaluate(variableDictionary_t vars) {

        iterator_type iter = expr_.begin();
        iterator_type end  = expr_.end();
        client::error_handler<iterator_type> error_handler(iter, end);
        client::code_gen::StackEvaluator     evaluator(error_handler);
        evaluator.registerVariables(vars);
        evaluator.registerFunctions(known_expr_funcs_);

        double result = 0.0;
        bool   valid  = false;
        if (evaluator(ast_)) {
            result = evaluator.result();
            valid  = true;
        }

        return boost::make_tuple(result, valid);
    }


private:

    typedef std::string::const_iterator iterator_type;
    client::ast::expression ast_;

    std::set<std::string> vars_;
    std::set<std::string> funcs_;

    std::string expr_;
    functionDictionary_t known_expr_funcs_;

    Expressions::OperatorType_t type_;

    void determineConstrOperator() {

        std::string op = expr_;

        if(boost::find_first(op, "=="))
            type_ = Expressions::EQ;
        else if(boost::find_first(op, "!="))
            type_ = Expressions::NOT_EQ;
        else if(boost::find_first(op, "<="))
            type_ = Expressions::INEQ_LHS_EQ;
        else if(boost::find_first(op, "<"))
            type_ = Expressions::INEQ_LHS;
        else if(boost::find_first(op, ">="))
            type_ = Expressions::INEQ_RHS_EQ;
        else if(boost::find_first(op,">"))
            type_ = Expressions::INEQ_RHS;
        else
            type_ = Expressions::NONE;
    }

    void parse() {

        iterator_type iter = expr_.begin();
        iterator_type end  = expr_.end();

        client::error_handler<iterator_type>       error_handler(iter, end);
        client::parser::expression<iterator_type>  expression(error_handler);
        client::parser::skipper<iterator_type>     skipper;
        client::code_gen::requirements             requirements(error_handler);

        bool success = phrase_parse(iter, end, expression, skipper, ast_);

        if (!success || iter != end) {
            std::cout << "Parsing failed!" << std::endl;
            std::string here = (iter != end ? std::string(iter, end): expr_);
            throw new OptPilotException("Expression::parse()",
                                         "Parsing failed here: " + here + "!");
        }

        // store the functions and variables required to evaluate this
        // expression
        if (requirements(ast_)) {
            vars_ = requirements.variables();
            funcs_ = requirements.functions();
        }
    }
};

#endif