//
// Struct StackEvaluator
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
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
#ifndef STACKEVALUATOR_HPP
#define STACKEVALUATOR_HPP

#include <cassert>
#include <map>
#include <string>

#include "ast.hpp"
#include "function.hpp"

namespace client { namespace code_gen
{
    struct StackEvaluator {

        typedef bool result_type;

        template <typename ErrorHandler>
        StackEvaluator(ErrorHandler&)
        {
        }

        double result() {
            assert(evaluation_stack_.size() == 1);
            client::function::argument_t res = evaluation_stack_.back();
            double result = std::get<double>(res);
            evaluation_stack_.pop_back();
            return result;
        }

        void registerFunction(std::string name,
                              client::function::type callback) {
            functions_.insert(client::function::named_t(name, callback));
        }

        void registerFunctions (
                std::map<std::string, client::function::type> functions) {
            functions_.insert(functions.begin(), functions.end());
        }

        void registerVariables(
                std::map<std::string, double> variableDictionary) {
            variableDictionary_ = variableDictionary;
        }

        // visitor
        bool operator()(ast::nil) { assert(false); return false; }
        bool operator()(unsigned int x);
        bool operator()(double x);
        bool operator()(bool x);
        bool operator()(ast::quoted_string const& x);
        bool operator()(ast::identifier const& x);
        bool operator()(ast::operation const& x);
        bool operator()(ast::unary const& x);
        bool operator()(ast::function_call const& x);
        bool operator()(ast::expression const& x);

    private:
        std::map<std::string, double> variableDictionary_;
        std::map<std::string, client::function::type> functions_;

        //our stack is conform to function call arguments
        client::function::arguments_t evaluation_stack_;
    };
}}

#endif