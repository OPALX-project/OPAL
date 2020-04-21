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
#ifndef COLUMN_HPP_
#define COLUMN_HPP_

#include "ast.hpp"
#include "skipper.hpp"
#include "error_handler.hpp"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>

#include <vector>

#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#define BOOST_SPIRIT_QI_DEBUG

namespace SDDS {
    struct column
    {
        enum attributes { NAME
                        , SYMBOL
                        , UNITS
                        , DESCRIPTION
                        , FORMAT_STRING
                        , TYPE
                        , FIELD_LENGTH
        };

        unsigned int order_m;
        boost::optional<std::string> name_m;
        boost::optional<std::string> units_m;
        boost::optional<std::string> description_m;
        boost::optional<ast::datatype> type_m;
        ast::columnData_t values_m;
        static unsigned int count_m;

        bool checkMandatories() const
        {
            return name_m && type_m;
        }

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply()
            {
                std::string attributeString;
                switch(A)
                {
                case SYMBOL:
                    attributeString = "symbol";
                    break;
                case FORMAT_STRING:
                    attributeString = "format_string";
                    break;
                case FIELD_LENGTH:
                    attributeString = "field_length";
                    break;
                default:
                    return true;
                }
                std::cerr << attributeString << " not supported yet" << std::endl;
                return false;
            }
        };

        template <typename Iterator, typename Skipper>
        bool parse(
            Iterator& first
          , Iterator last
          , Skipper const& skipper)
        {
            switch(*this->type_m) {
            case ast::FLOAT:
            {
                float f = 0.0;
                boost::spirit::qi::float_type float_;
                if (phrase_parse(first, last, float_, skipper, f)) {
                    this->values_m.push_back(f);
                    return true;
                }
                break;
            }
            case ast::DOUBLE:
            {
                double d = 0.0;
                boost::spirit::qi::double_type double_;
                if (phrase_parse(first, last, double_, skipper, d)) {
                    this->values_m.push_back(d);
                    return true;
                }
                break;
            }
            case ast::SHORT:
            {
                short s = 0;
                boost::spirit::qi::short_type short_;
                if (phrase_parse(first, last, short_, skipper, s)) {
                    this->values_m.push_back(s);
                    return true;
                }
                break;
            }
            case ast::LONG:
            {
                long l = 0;
                boost::spirit::qi::long_type long_;
                if (phrase_parse(first, last, long_, skipper, l)) {
                    this->values_m.push_back(l);
                    return true;
                }
                break;
            }
            case ast::CHARACTER:
            {
                char c = '\0';
                boost::spirit::qi::char_type char_;
                if (phrase_parse(first, last, char_, skipper, c)) {
                    this->values_m.push_back(c);
                    return true;
                }
                break;
            }
            case ast::STRING:
            {
                std::string s("");
                parser::qstring<Iterator, Skipper> qstring;
                if (phrase_parse(first, last, qstring, skipper, s)) {
                    this->values_m.push_back(s);
                    return true;
                }
                break;
            }
            }
            return false;
        }
    };

    struct columnList: std::vector<column> {};

    template <typename Iterator>
    struct columnOrder
    {
        template <typename, typename>
        struct result { typedef void type; };

        void operator()(column& col, Iterator) const
        {
            col.order_m = column::count_m ++;
        }
    };

    inline std::ostream& operator<<(std::ostream& out, const column& col) {
        if (col.name_m) out << "name = " << *col.name_m << ", ";
        if (col.type_m) out << "type = " << *col.type_m << ", ";
        if (col.units_m) out << "units = " << *col.units_m << ", ";
        if (col.description_m) out << "description = " << *col.description_m << ", ";
        out << "order = " << col.order_m;

        return out;
    }
}

BOOST_FUSION_ADAPT_STRUCT(
    SDDS::column,
    (boost::optional<std::string>, name_m)
    (boost::optional<SDDS::ast::datatype>, type_m)
    (boost::optional<std::string>, units_m)
    (boost::optional<std::string>, description_m)
    (SDDS::ast::variant_t, value_m)
)

namespace SDDS { namespace parser
{
    namespace qi = boost::spirit::qi;
    namespace ascii = boost::spirit::ascii;
    namespace phx = boost::phoenix;

    ///////////////////////////////////////////////////////////////////////////////
    //  The expression grammar
    ///////////////////////////////////////////////////////////////////////////////
    template <typename Iterator>
    struct column_parser: qi::grammar<Iterator, column(), skipper<Iterator> >
    {
        column_parser(error_handler<Iterator> & _error_handler);

        qi::rule<Iterator, std::string(), skipper<Iterator> > string, quoted_string, units;
        qi::rule<Iterator, std::string(), skipper<Iterator> > column_name, column_units,
                column_description, column_symbol, column_format;
        qi::rule<Iterator, ast::datatype(), skipper<Iterator> > column_type;
        qi::rule<Iterator, column(), skipper<Iterator> > start;
        qi::rule<Iterator, long(), skipper<Iterator> > column_field;
        qi::rule<Iterator, ast::nil(), skipper<Iterator> > column_unsupported_pre,
                column_unsupported_post;
        qi::symbols<char, ast::datatype> datatype;
    };
}}
#endif /* COLUMN_HPP_ */