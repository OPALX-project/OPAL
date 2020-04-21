//
// Struct file
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
#ifndef FILE_HPP_
#define FILE_HPP_

#include "ast.hpp"
#include "skipper.hpp"
#include "error_handler.hpp"
#include "version.hpp"
#include "description.hpp"
#include "parameter.hpp"
#include "column.hpp"
#include "data.hpp"
#include "associate.hpp"
#include "array.hpp"
#include "include.hpp"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>

#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#define BOOST_SPIRIT_QI_DEBUG

namespace SDDS {
    struct file
    {
        version sddsVersion_m;                              // 0
        boost::optional<description> sddsDescription_m;     // 1
        parameterList sddsParameters_m;                     // 2
        columnList sddsColumns_m;                           // 3
        data sddsData_m;                                    // 4
        associateList sddsAssociates_m;                     // 5
        arrayList sddsArrays_m;                             // 6
        includeList sddsIncludes_m;                         // 7

        void clear() {
            sddsParameters_m.clear();
            sddsColumns_m.clear();
            sddsAssociates_m.clear();
            sddsArrays_m.clear();
            sddsIncludes_m.clear();
        }
    };

    inline std::ostream& operator<<(std::ostream& out, const file& fh) {
        out << fh.sddsVersion_m << "\n";
        if (fh.sddsDescription_m) out << *fh.sddsDescription_m << "\n";
        for (const parameter& param: fh.sddsParameters_m) {
            out << param << "\n";
        }
        for (const column& col: fh.sddsColumns_m) {
            out << col << "\n";
        }
        out << fh.sddsData_m << "\n";

//        if (fh.sddsAssociates_m) out << *fh.sddsAssociates_m << "\n";
//        if (fh.sddsArrays_m) out << *fh.sddsArrays_m << "\n";
//        if (fh.sddsIncludes_m) out << *fh.sddsIncludes_m << "\n";

        out << std::endl;

        return out;
    }
}

BOOST_FUSION_ADAPT_STRUCT(
    SDDS::file,
    (SDDS::version, sddsVersion_m)
    (boost::optional<SDDS::description>, sddsDescription_m)
    (SDDS::parameterList, sddsParameters_m)
    (SDDS::columnList, sddsColumns_m)
    (SDDS::data, sddsData_m)
    (SDDS::associateList, sddsAssociates_m)
    (SDDS::arrayList, sddsArrays_m)
    (SDDS::includeList, sddsIncludes_m)
)

namespace SDDS { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phx = boost::phoenix;

    template <typename Iterator>
    struct file_parser: qi::grammar<Iterator, file(), skipper<Iterator> >
    {
        file_parser(error_handler<Iterator> & _error_handler);

        version_parser<Iterator> version_m;
        description_parser<Iterator> description_m;
        parameter_parser<Iterator> parameter_m;
        column_parser<Iterator> column_m;
        data_parser<Iterator> data_m;
        associate_parser<Iterator> associate_m;
        array_parser<Iterator> array_m;
        include_parser<Iterator> include_m;
        qi::rule<Iterator, file(), skipper<Iterator> > start;
    };
}}
#endif /* FILE_HPP_ */