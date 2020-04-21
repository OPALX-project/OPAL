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
#ifndef ARRAY_DEF_HPP_
#define ARRAY_DEF_HPP_

#include "array.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    array_parser<Iterator>::array_parser(error_handler<Iterator> & _error_handler):
        array_parser::base_type(start)
    {
        using qi::on_error;
        using qi::fail;
        using phx::function;
        typedef function<error_handler<Iterator> > error_handler_function;

        // qi::_1_type _1;
        qi::_3_type _3;
        qi::_4_type _4;
        qi::lexeme_type lexeme;
        qi::char_type char_;
        qi::lit_type lit;
        qi::long_type long_;
        qi::short_type short_;
        qi::alpha_type alpha;
        qi::alnum_type alnum;
        // qi::_val_type _val;
        qi::_pass_type _pass;
        qi::eps_type eps;

        quoted_string %= lexeme['"' >> +(char_ - '"') >> '"'];
        string %= quoted_string
                | lexeme[(alpha | char_("@:#+-%._$&/")) >> *(alnum | char_("@:#+-%._$&/"))];
        units %= lexeme[alpha >> *(alpha | '/')]
                 | lit('1');

         arraytype.add
             ("float", ast::FLOAT)
             ("double", ast::DOUBLE)
             ("short", ast::SHORT)
             ("long", ast::LONG)
             ("character", ast::CHARACTER)
             ("string", ast::STRING)
             ;

        array_name = lit("mode") > '=' > string;
        array_symbol = lit("symbol") > '=' > string;
        array_units = lit("units") > '=' > units;
        array_description = lit("description") > '=' > string;
        array_format = lit("format_string") > '=' > string;
        array_group = lit("group_name") > '=' > string;
        array_type = lit("type") > '=' > arraytype;
        array_field = lit("field_length") > '=' > long_;
        array_dimensions = lit("dimensions") > '=' > long_;

        auto complainArray = phx::bind(&array::complainUnsupported<array::ARRAY>::apply);
        auto complainName = phx::bind(&array::complainUnsupported<array::NAME>::apply);
        auto complainSymbol= phx::bind(&array::complainUnsupported<array::SYMBOL>::apply);
        auto complainUnits = phx::bind(&array::complainUnsupported<array::UNITS>::apply);
        auto complainDescription = phx::bind(&array::complainUnsupported<array::DESCRIPTION>::apply);
        auto complainFormat= phx::bind(&array::complainUnsupported<array::FORMAT_STRING>::apply);
        auto complainGroup = phx::bind(&array::complainUnsupported<array::GROUP_NAME>::apply);
        auto complainType = phx::bind(&array::complainUnsupported<array::TYPE>::apply);
        auto complainField = phx::bind(&array::complainUnsupported<array::FIELD_LENGTH>::apply);
        auto complainDimensions = phx::bind(&array::complainUnsupported<array::DIMENSIONS>::apply);

        array_unsupported_pre =
                ((',' > array_name[_pass = complainName])
                 ^ (',' > array_symbol[_pass = complainSymbol])
                 ^ (',' > array_units[_pass = complainUnits])
                 ^ (',' > array_description[_pass = complainDescription])
                 ^ (',' > array_format[_pass = complainFormat])
                 ^ (',' > array_group[_pass = complainGroup])
                 ^ (',' > array_type[_pass = complainType])
                 ^ (',' > array_field[_pass = complainField])
                 ^ (',' > array_dimensions[_pass = complainDimensions])
        );
        array_unsupported_post =
                ((array_name[_pass = complainName] > ',')
                 ^ (array_symbol[_pass = complainSymbol] > ',')
                 ^ (array_units[_pass = complainUnits] > ',')
                 ^ (array_description[_pass = complainDescription] > ',')
                 ^ (array_format[_pass = complainFormat] > ',')
                 ^ (array_group[_pass = complainGroup] > ',')
                 ^ (array_type[_pass = complainType] > ',')
                 ^ (array_field[_pass = complainField] > ',')
                 ^ (array_dimensions[_pass = complainDimensions] > ',')
        );

        start =
                lit("&array")[_pass = complainArray]
                > lit("&end");

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                  std::string("Error! Expecting "), _4, _3));
    }
}}
#endif /* ARRAY_DEF_HPP_ */