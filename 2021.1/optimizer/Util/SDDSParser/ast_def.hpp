//
// Namespace ast
//
// Copyright (c) 2015, Christof Metzger-Kraus, Helmholtz-Zentrum Berlin
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
#ifndef AST_DEF_HPP_
#define AST_DEF_HPP_

#include "ast.hpp"

namespace SDDS { namespace parser
{

    template <typename Iterator, typename Skipper>
    string<Iterator, Skipper>::string():
        string::base_type(start)
    {
        qi::char_type char_;
        qi::eol_type eol;
        qi::lexeme_type lexeme;

        start %= lexeme[+(char_-eol)];
    }


    template <typename Iterator, typename Skipper>
    qstring<Iterator, Skipper>::qstring():
        qstring::base_type(start)
    {
        qi::char_type char_;
        qi::eol_type eol;
        qi::lexeme_type lexeme;

        start %= lexeme['"' >> +(char_ - (eol|'"')) >> '"']
               | lexeme[+(char_ - (eol|' '))];
    }
}}

#endif // AST_DEF_HPP_
