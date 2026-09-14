//
// Namespace requirements
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
#ifndef REQUIREMENTS_HPP
#define REQUIREMENTS_HPP

#include "ast.hpp"

#include <cassert>
#include <set>

namespace client { namespace code_gen
{
    struct requirements
    {
        typedef bool result_type;

        template <typename ErrorHandler>
        requirements(ErrorHandler&)
        {
        }

        bool operator()(ast::nil) { assert(false); return false; }
        bool operator()(unsigned int /*x*/) { return true; }
        bool operator()(double /*x*/) { return true; }
        bool operator()(bool /*x*/) { return true; }
        bool operator()(ast::quoted_string const & /*x*/) { return true; }

        bool operator()(ast::operation const& x) {
            if (!ast::apply_visitor(*this, x.operand_))
                return false;
            return true;
        }

        bool operator()(ast::unary const& x) {
            if (!ast::apply_visitor(*this, x.operand_))
                return false;
            return true;
        }

        bool operator()(ast::identifier const& x) {
            variables_.insert(x.name);

            return true;
        }

        bool operator()(ast::function_call const& x) {
            functions_.insert(x.function_name.name);

            for(ast::function_call_argument const& arg: x.args) {
                if (!ast::apply_visitor(*this, arg))
                    return false;
                //if (!(*this)(arg))
                    //return false;
            }
            return true;
        }

        bool operator()(ast::expression const& x) {

            if (!ast::apply_visitor(*this, x.first))
                return false;

            for (ast::operation const& oper: x.rest) {
                if (!(*this)(oper))
                    return false;
        }

        return true;
    }

    std::set<std::string> variables() { return variables_; }
    std::set<std::string> functions() { return functions_; }

    private:
        std::set<std::string> variables_;
        std::set<std::string> functions_;
    };
}}

#endif