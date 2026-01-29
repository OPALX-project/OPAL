//
// Struct skipper
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
#ifndef SKIPPER_HPP_
#define SKIPPER_HPP_

#include <boost/spirit/include/qi.hpp>

namespace SDDS { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;

    ///////////////////////////////////////////////////////////////////////////////
    //  The skipper grammar
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct skipper : qi::grammar<Iterator>
    {
        skipper() : skipper::base_type(start)
        {
        	qi::eol_type eol;
        	qi::eoi_type eoi;
        	qi::char_type char_;
            ascii::space_type space;

            start =
                    space
				|   "!" >> *(char_ - eol) >> (eol|eoi)       // comments
                ;
        }

        qi::rule<Iterator> start;
    };

    template <typename Iterator>
    struct listskipper : qi::grammar<Iterator>
    {
        listskipper() : listskipper::base_type(start)
        {
            qi::char_type char_;
            ascii::space_type space;

            start =
                    space
                |   char_(',')
                ;
        }

        qi::rule<Iterator> start;
    };

}}

#endif


