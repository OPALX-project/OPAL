//
// Struct column
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
#ifndef COLUMN_DEF_HPP_
#define COLUMN_DEF_HPP_

#include "column.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    column_parser<Iterator>::column_parser(error_handler<Iterator> & _error_handler):
        column_parser::base_type(start)
    {
        using qi::on_error;
        using qi::on_success;
        using qi::fail;
        using phx::function;
        typedef function<error_handler<Iterator> > error_handler_function;

        // qi::_0_type _0;
        qi::_1_type _1;
        qi::_3_type _3;
        qi::_4_type _4;
        qi::lexeme_type lexeme;
        qi::char_type char_;
        qi::long_type long_;
        qi::lit_type lit;
        qi::alpha_type alpha;
        qi::alnum_type alnum;
        qi::_val_type _val;
        qi::_pass_type _pass;
        qi::eps_type eps;

        quoted_string %= lexeme['"' >> +(char_ - '"') >> '"'];
        string %= quoted_string
                | lexeme[(alpha | char_("@:#+-%._$&/")) >> *(alnum | char_("@:#+-%._$&/"))];
        units %= char_('1') | lexeme[alpha >> *(alpha | char_('/'))];

        datatype.add
            ("float", ast::FLOAT)
            ("double", ast::DOUBLE)
            ("short", ast::SHORT)
            ("long", ast::LONG)
            ("character", ast::CHARACTER)
            ("string", ast::STRING)
            ;

        column_name = lit("name") > '=' > string;
        column_units %= lit("units") > '=' > units;
        column_description %= lit("description") > '=' > string;
        column_type %= lit("type") > '=' > datatype;
        column_symbol = lit("symbol") > '=' > string;
        column_format = lit("format_string") > '=' > string;
        column_field = lit("field_length") > '=' > long_;

        auto complainSymbol = phx::bind(&column::complainUnsupported<column::SYMBOL>::apply);
        auto complainFormat = phx::bind(&column::complainUnsupported<column::FORMAT_STRING>::apply);
        auto complainField = phx::bind(&column::complainUnsupported<column::FIELD_LENGTH>::apply);

        column_unsupported_pre =
                ((',' > column_symbol[_pass = complainSymbol])
                 ^ (',' > column_format[_pass = complainFormat])
                 ^ (',' > column_field[_pass = complainField])
                 );
        column_unsupported_post =
                ((column_symbol[_pass = complainSymbol] > ',')
                 ^ (column_format[_pass = complainFormat] > ',')
                 ^ (column_field[_pass = complainField] > ',')
                 );

        start =
                lit("&column")
                > -column_unsupported_post
                >> ((column_name[phx::at_c<0>(_val) = _1]
                     >> ((',' >> column_type[phx::at_c<1>(_val) = _1])
                        ^(',' >> column_units[phx::at_c<2>(_val) = _1])
                        ^(',' >> column_description[phx::at_c<3>(_val) = _1])
                        ))
                    |(column_type[phx::at_c<1>(_val) = _1]
                      >> ((',' >> column_name[phx::at_c<0>(_val) = _1])
                         ^(',' >> column_units[phx::at_c<2>(_val) = _1])
                         ^(',' >> column_description[phx::at_c<3>(_val) = _1])
                         ))
                    |(column_units[phx::at_c<2>(_val) = _1]
                      >> ((',' >> column_type[phx::at_c<1>(_val) = _1])
                         ^(',' >> column_name[phx::at_c<0>(_val) = _1])
                         ^(',' >> column_description[phx::at_c<3>(_val) = _1])
                         ))
                    |(column_description[phx::at_c<3>(_val) = _1]
                      >> ((',' >> column_type[phx::at_c<1>(_val) = _1])
                         ^(',' >> column_units[phx::at_c<2>(_val) = _1])
                         ^(',' >> column_name[phx::at_c<0>(_val) = _1])
                         ))
                   )
                > -column_unsupported_pre
                > lit("&end")
                >> eps[_pass = phx::bind(&column::checkMandatories, _val)];

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                  std::string("Error! Expecting "), _4, _3));

    }
}}
#endif /* COLUMN_DEF_HPP_ */