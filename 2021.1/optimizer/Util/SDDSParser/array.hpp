//
// Struct array
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
#ifndef ARRAY_HPP_
#define ARRAY_HPP_

#include "ast.hpp"
#include "skipper.hpp"
#include "error_handler.hpp"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>

#include <list>

#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#define BOOST_SPIRIT_QI_DEBUG

namespace SDDS {
    struct array
    {
        enum attributes { NAME
                         , SYMBOL
                         , UNITS
                         , DESCRIPTION
                         , FORMAT_STRING
                         , GROUP_NAME
                         , TYPE
                         , FIELD_LENGTH
                         , DIMENSIONS
                         , ARRAY
         };

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply()
            {
                std::string attributeString;
                switch(A)
                {
                case NAME:
                    attributeString = "name";
                    break;
                case SYMBOL:
                    attributeString = "symbol";
                    break;
                case UNITS:
                    attributeString = "units";
                    break;
                case DESCRIPTION:
                    attributeString = "description";
                    break;
                case FORMAT_STRING:
                    attributeString = "format_string";
                    break;
                case GROUP_NAME:
                    attributeString = "group_name";
                    break;
                case TYPE:
                    attributeString = "type";
                    break;
                case FIELD_LENGTH:
                    attributeString = "field_length";
                    break;
                case DIMENSIONS:
                    attributeString = "dimensions";
                    break;
                case ARRAY:
                    attributeString = "array";
                    break;
                default:
                    return true;
                }
                std::cerr << attributeString << " not supported yet" << std::endl;
                return false;
            }
        };
    };

    struct arrayList: std::list<array> {};

    inline std::ostream& operator<<(std::ostream& out, const array& ) {
        return out;
    }
}

namespace SDDS { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phx = boost::phoenix;

    template <typename Iterator>
    struct array_parser: qi::grammar<Iterator, array(), skipper<Iterator> >
    {
        array_parser(error_handler<Iterator> & _error_handler);

        qi::rule<Iterator, array(), skipper<Iterator> > start;
        qi::rule<Iterator, std::string(), skipper<Iterator> > array_name,
                array_symbol, array_units, array_description, array_format,
                array_group, units, string, quoted_string;
        qi::rule<Iterator, long(), skipper<Iterator> > array_field,
                array_dimensions;
        qi::rule<Iterator, ast::datatype, skipper<Iterator> > array_type;
        qi::rule<Iterator, ast::nil(), skipper<Iterator> > array_unsupported_pre,
                array_unsupported_post;
        qi::symbols<char, ast::datatype> arraytype;
    };
}}
#endif /* ARRAY_HPP_ */
