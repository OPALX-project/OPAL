//
// Struct parameter
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
#ifndef PARAMETER_DEF_HPP_
#define PARAMETER_DEF_HPP_

#include "parameter.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    parameter_parser<Iterator>::parameter_parser(error_handler<Iterator> & _error_handler):
        parameter_parser::base_type(start)
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
        qi::lit_type lit;
        qi::long_type long_;
        qi::alpha_type alpha;
        qi::alnum_type alnum;
        qi::_val_type _val;
        qi::_pass_type _pass;
        qi::eps_type eps;

        quoted_string %= lexeme['"' >> +(char_ - '"') >> '"'];
        string %= quoted_string
                | lexeme[(alpha | char_("@:#+-%._$&/")) >> *(alnum | char_("@:#+-%._$&/"))];
        units %= lexeme[alpha >> *(alpha | char_('/'))]
                | char_('1');

        datatype.add
            ("float", ast::FLOAT)
            ("double", ast::DOUBLE)
            ("short", ast::SHORT)
            ("long", ast::LONG)
            ("character", ast::CHARACTER)
            ("string", ast::STRING)
            ;

        parameter_name = lit("name") > '=' > string;
        parameter_units %= lit("units") > '=' > units;
        parameter_description %= lit("description") > '=' > string;
        parameter_type %= lit("type") > '=' > datatype;
        parameter_symbol = lit("symbol") > '=' > string;
        parameter_format = lit("format_string") > '=' > string;
        parameter_fixed = lit("fixed_value") > long_;

        auto complainSymbol = phx::bind(&parameter::complainUnsupported<parameter::SYMBOL>::apply);
        auto complainFormat = phx::bind(&parameter::complainUnsupported<parameter::FORMAT_STRING>::apply);
        auto complainFixed = phx::bind(&parameter::complainUnsupported<parameter::FIXED_VALUE>::apply);

        parameter_unsupported_pre =
                ((',' > parameter_symbol[_pass = complainSymbol])
               ^ (',' > parameter_format[_pass = complainFormat])
               ^ (',' > parameter_fixed[_pass = complainFixed])
                 );
        parameter_unsupported_post =
                ((parameter_symbol[_pass = complainSymbol] > ',')
               ^ (parameter_format[_pass = complainFormat] > ',')
               ^ (parameter_fixed[_pass = complainFixed] > ',')
                 );

        start =
                lit("&parameter")
                > -parameter_unsupported_post
                >> ((parameter_name[phx::at_c<0>(_val) = _1]
                     >> ((',' >> parameter_type[phx::at_c<1>(_val) = _1])
                        ^(',' >> parameter_units[phx::at_c<2>(_val) = _1])
                        ^(',' >> parameter_description[phx::at_c<3>(_val) = _1])
                        ))
                    |(parameter_type[phx::at_c<1>(_val) = _1]
                      >> ((',' >> parameter_name[phx::at_c<0>(_val) = _1])
                         ^(',' >> parameter_units[phx::at_c<2>(_val) = _1])
                         ^(',' >> parameter_description[phx::at_c<3>(_val) = _1])
                         ))
                    |(parameter_units[phx::at_c<2>(_val) = _1]
                      >> ((',' >> parameter_type[phx::at_c<1>(_val) = _1])
                         ^(',' >> parameter_name[phx::at_c<0>(_val) = _1])
                         ^(',' >> parameter_description[phx::at_c<3>(_val) = _1])
                         ))
                    |(parameter_description[phx::at_c<3>(_val) = _1]
                      >> ((',' >> parameter_type[phx::at_c<1>(_val) = _1])
                         ^(',' >> parameter_units[phx::at_c<2>(_val) = _1])
                         ^(',' >> parameter_name[phx::at_c<0>(_val) = _1])
                         ))
                   )
                > -parameter_unsupported_pre
                > lit("&end")
                >> eps[_pass = phx::bind(&parameter::checkMandatories, _val)];

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                  std::string("Error! Expecting "), _4, _3));

    }
}}
#endif /* PARAMETER_DEF_HPP_ */