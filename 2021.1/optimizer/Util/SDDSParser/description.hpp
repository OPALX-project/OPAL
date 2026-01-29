//
// Struct description
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
#ifndef DESCRIPTION_HPP_
#define DESCRIPTION_HPP_

#include "ast.hpp"
#include "skipper.hpp"
#include "error_handler.hpp"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>

#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#define BOOST_SPIRIT_QI_DEBUG

namespace SDDS {
    struct description
    {
        boost::optional<std::string> text_m;
        boost::optional<std::string> content_m;

    };

    inline std::ostream& operator<<(std::ostream& out, const description& desc) {
        if (desc.text_m) out << "text = " << *desc.text_m << ", ";
        if (desc.content_m) out << "content = " << *desc.content_m;
        return out;
    }
}

BOOST_FUSION_ADAPT_STRUCT(
    SDDS::description,
    (boost::optional<std::string>, text_m)
    (boost::optional<std::string>, content_m)
)

namespace SDDS { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phx = boost::phoenix;

    template <typename Iterator>
    struct description_parser: qi::grammar<Iterator, description(), skipper<Iterator> >
    {
        description_parser(error_handler<Iterator> & _error_handler);

        qi::rule<Iterator, description(), skipper<Iterator> > start;
        qi::rule<Iterator, std::string(), skipper<Iterator> > quoted_string,
                description_text, description_content;
    };
}}
#endif /* DESCRIPTION_HPP_ */
