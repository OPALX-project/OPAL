//
// Namespace parser expression
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
#ifndef EXPRESSION_HPP
#define EXPRESSION_HPP

#include "annotation.hpp"
#include "ast.hpp"
#include "error_handler.hpp"

#include <list>
#include <string>

namespace client { namespace parser
{
    ///////////////////////////////////////////////////////////////////////////////
    //  The expression grammar
    //
    //  This is a hand-written recursive-descent parser (replacing the former
    //  Boost.Spirit Qi grammar) that builds the same ast::expression tree.
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct expression
    {
        explicit expression(error_handler<Iterator>& error_handler);

        /// Parses an expression starting at `first`. On success, `first` is
        /// advanced past the consumed input and `attr` is filled in; on
        /// failure `first` is left unchanged and `false` is returned.
        bool parse(Iterator& first, Iterator last, ast::expression& attr);

    private:
        error_handler<Iterator>&  error_handler_;
        annotation<Iterator>      annotation_;

        void skip_ws(Iterator& it, Iterator last) const;
        bool match(Iterator& it, Iterator last, char c) const;
        bool match_word(Iterator& it, Iterator last, char const* word) const;
        bool match_op_token(Iterator& it, Iterator last, char const* token) const;
        bool try_identifier_token(Iterator& it, Iterator last, std::string& word) const;

        template <typename NextLevel, typename OpMatcher>
        bool parse_level(Iterator& it, Iterator last, ast::expression& out,
                          NextLevel next, OpMatcher match_op);

        bool parse_expr(Iterator& it, Iterator last, ast::expression& out);
        bool parse_logical_or(Iterator& it, Iterator last, ast::expression& out);
        bool parse_logical_and(Iterator& it, Iterator last, ast::expression& out);
        bool parse_equality(Iterator& it, Iterator last, ast::expression& out);
        bool parse_relational(Iterator& it, Iterator last, ast::expression& out);
        bool parse_additive(Iterator& it, Iterator last, ast::expression& out);
        bool parse_multiplicative(Iterator& it, Iterator last, ast::expression& out);

        bool parse_unary(Iterator& it, Iterator last, ast::operand& out);
        bool parse_primary(Iterator& it, Iterator last, ast::operand& out);
        bool parse_constant(Iterator& it, Iterator last, ast::operand& out);

        bool parse_argument_list(Iterator& it, Iterator last, std::list<ast::function_call_argument>& out);
        bool parse_call_argument(Iterator& it, Iterator last, ast::function_call_argument& out);
        bool parse_quoted_string(Iterator& it, Iterator last, std::string& out);

        static ast::operand to_operand(ast::expression&& e);
    };
}}

#endif