//
// Namespace annotation
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
#ifndef ANNOTATION_HPP
#define ANNOTATION_HPP

#include <vector>

#include "ast.hpp"

namespace client
{
    ///////////////////////////////////////////////////////////////////////////////
    //  The annotation handler links the AST to a map of iterator positions
    //  for the purpose of subsequent semantic error handling when the
    //  program is being compiled.
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct annotation
    {
        template <typename>
        struct result { typedef void type; };

        std::vector<Iterator>& iters;
        annotation(std::vector<Iterator>& iters)
          : iters(iters) {}

        struct set_id
        {
            typedef void result_type;

            int id;
            set_id(int id) : id(id) {}

            void operator()(ast::function_call& x) const
            {
                x.function_name.id = id;
            }

            void operator()(ast::identifier& x) const
            {
                x.id = id;
            }

            template <typename T>
            void operator()(T& /*x*/) const
            {
                // no-op
            }
        };

        void operator()(ast::operand& ast, Iterator pos) const
        {
            int id = iters.size();
            iters.push_back(pos);
            ast::apply_visitor(set_id(id), ast);
        }

        void operator()(ast::identifier& ast, Iterator pos) const
        {
            int id = iters.size();
            iters.push_back(pos);
            ast.id = id;
        }
    };
}

#endif