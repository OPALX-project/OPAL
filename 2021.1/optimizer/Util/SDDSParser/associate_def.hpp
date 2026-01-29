//
// Struct associate
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
#ifndef ASSOCIATE_DEF_HPP_
#define ASSOCIATE_DEF_HPP_

#include "associate.hpp"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_bind.hpp>

namespace SDDS { namespace parser
{
    template <typename Iterator>
    associate_parser<Iterator>::associate_parser(error_handler<Iterator> & _error_handler):
        associate_parser::base_type(start)
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
                | lexeme[(alpha | '_') >> *(alnum | '_')];

        associate_name = lit("mode") > '=' > string;
        associate_filename = lit("filename") > '=' > string;
        associate_path = lit("path") > '=' > string;
        associate_description = lit("description") > '=' > string;
        associate_contents= lit("contents") > '=' > string;
        associate_sdds= lit("column_major_order") > '=' > long_;

        auto complainAssociate = phx::bind(&associate::complainUnsupported<associate::ASSOCIATE>::apply);
        auto complainName = phx::bind(&associate::complainUnsupported<associate::NAME>::apply);
        auto complainFilename = phx::bind(&associate::complainUnsupported<associate::FILENAME>::apply);
        auto complainPath = phx::bind(&associate::complainUnsupported<associate::PATH>::apply);
        auto complainDescription = phx::bind(&associate::complainUnsupported<associate::DESCRIPTION>::apply);
        auto complainContents = phx::bind(&associate::complainUnsupported<associate::CONTENTS>::apply);
        auto complainSDDS = phx::bind(&associate::complainUnsupported<associate::SDDS>::apply);

        associate_unsupported_pre =
                ((',' > associate_name[_pass = complainName])
                 ^ (',' > associate_filename[_pass = complainFilename])
                 ^ (',' > associate_path[_pass = complainPath])
                 ^ (',' > associate_description[_pass = complainDescription])
                 ^ (',' > associate_contents[_pass = complainContents])
                 ^ (',' > associate_sdds[_pass = complainSDDS]));
        associate_unsupported_post =
                ((associate_name[_pass = complainName] > ',')
                 ^ (associate_filename[_pass = complainFilename] > ',')
                 ^ (associate_path[_pass = complainPath] > ',')
                 ^ (associate_description[_pass = complainDescription] > ',')
                 ^ (associate_contents[_pass = complainContents] > ',')
                 ^ (associate_sdds[_pass = complainSDDS] > ','));

        start =
                lit("&associate")[_pass = complainAssociate]
                > lit("&end");

        BOOST_SPIRIT_DEBUG_NODES(
            (start)
        )

        on_error<fail>(start,
            error_handler_function(_error_handler)(
                  std::string("Error! Expecting "), _4, _3));
    }
}}
#endif /* ASSOCIATE_HPP_ */