//
// Struct data
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
#ifndef DATA_DEF_HPP_
#define DATA_DEF_HPP_

#include "data.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    data_parser<Iterator>::data_parser(error_handler<Iterator> & _error_handler):
        data_parser::base_type(start)
    {
        using qi::on_error;
        using qi::fail;
        using phx::function;
        typedef function<error_handler<Iterator> > error_handler_function;

        qi::_1_type _1;
        qi::_3_type _3;
        qi::_4_type _4;
        qi::lit_type lit;
        qi::long_type long_;
        qi::short_type short_;
        qi::_val_type _val;
        qi::_pass_type _pass;
        qi::eps_type eps;

        datamode.add
            ("ascii", ast::ASCII)
            ("binary", ast::BINARY)
            ;
        dataendian.add
            ("big", ast::BIGENDIAN)
            ("little", ast::LITTLEENDIAN)
            ;

        data_mode = lit("mode") > '=' > datamode;
        data_lines = lit("lines_per_row") > '=' > long_;
        data_row = lit("no_row_counts") > '=' > long_;
        data_fixed = lit("fixed_row_count") > '=' > long_;
        data_additional = lit("additional_header_lines") > '=' > long_;
        data_column = lit("column_major_order") > '=' > short_;
        data_endian = lit("endian") > '=' > dataendian;

        auto complainLines = phx::bind(&data::complainUnsupported<data::LINES_PER_ROW>::apply);
        auto complainFixed = phx::bind(&data::complainUnsupported<data::FIXED_ROW_COUNT>::apply);
        auto complainAdditional = phx::bind(&data::complainUnsupported<data::ADDITIONAL_HEADER_LINES>::apply);
        auto complainColumn = phx::bind(&data::complainUnsupported<data::COLUMN_MAJOR_ORDER>::apply);
        auto complainEndian = phx::bind(&data::complainUnsupported<data::ENDIAN>::apply);

        data_unsupported_pre =
                ((',' > data_lines[_pass = complainLines])
                 ^ (',' > data_fixed[_pass = complainFixed])
                 ^ (',' > data_additional[_pass = complainAdditional])
                 ^ (',' > data_column[_pass = complainColumn])
                 ^ (',' > data_endian[_pass = complainEndian]));
        data_unsupported_post =
                ((data_lines[_pass = complainLines] > ',')
                 ^ (data_fixed[_pass = complainFixed] > ',')
                 ^ (data_additional[_pass = complainAdditional] > ',')
                 ^ (data_column[_pass = complainColumn] > ',')
                 ^ (data_endian[_pass = complainEndian] > ','));

        start =
                lit("&data")
                > -data_unsupported_post
                > ((data_mode[phx::at_c<0>(_val) = _1] >> ',' >> data_row[phx::at_c<1>(_val) = _1])
                 | (data_row[phx::at_c<1>(_val) = _1] >> ',' >> data_mode[phx::at_c<0>(_val) = _1]))
                > -data_unsupported_pre
                > lit("&end")
                >> eps[_pass = phx::bind(&data::isASCII, _val)];

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                                                   std::string("Error! Expecting "), _4, _3));
    }
}}
#endif /* DATA_DEF_HPP_ */