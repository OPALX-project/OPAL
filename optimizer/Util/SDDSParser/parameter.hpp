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
#ifndef PARAMETER_HPP_
#define PARAMETER_HPP_

#include "ast.hpp"
#include "skipper.hpp"
#include "error_handler.hpp"

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/optional.hpp>

#include <list>

#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#define BOOST_SPIRIT_QI_DEBUG

namespace SDDS {
    struct parameter
    {
        enum attributes { NAME
                        , SYMBOL
                        , UNITS
                        , DESCRIPTION
                        , FORMAT_STRING
                        , TYPE
                        , FIXED_VALUE
         };

        unsigned int order_m;
        boost::optional<std::string> name_m;
        boost::optional<std::string> units_m;
        boost::optional<std::string> description_m;
        boost::optional<ast::datatype> type_m;
        ast::variant_t value_m;
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
                case FIXED_VALUE:
                    attributeString = "fixed_value";
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
                boost::spirit::qi::float_type float_;
                return phrase_parse(first, last, float_, skipper, this->value_m);
            }
            case ast::DOUBLE:
            {
                boost::spirit::qi::double_type double_;
                return phrase_parse(first, last, double_, skipper, this->value_m);
            }
            case ast::SHORT:
            {
                boost::spirit::qi::short_type short_;
                return phrase_parse(first, last, short_, skipper, this->value_m);
            }
            case ast::LONG:
            {
                boost::spirit::qi::long_type long_;
                return phrase_parse(first, last, long_, skipper, this->value_m);
            }
            case ast::CHARACTER:
            {
                boost::spirit::qi::char_type char_;
                return phrase_parse(first, last, char_, skipper, this->value_m);
            }
            case ast::STRING:
            {
                parser::string<Iterator, Skipper> string;
                return phrase_parse(first, last, string, skipper, this->value_m);
            }
            }
            return false;
        }
    };

    struct parameterList : std::vector<parameter> {};

    template <typename Iterator>
    struct parameterOrder
    {
        template <typename, typename>
        struct result { typedef void type; };

        void operator()(parameter& param, Iterator) const
        {
            param.order_m = parameter::count_m ++;
        }
    };

    inline std::ostream& operator<<(std::ostream& out, const parameter& param) {
        if (param.name_m) out << "name = " << *param.name_m << ", ";
        if (param.type_m) out << "type = " << *param.type_m << ", ";
        if (param.units_m) out << "units = " << *param.units_m << ", ";
        if (param.description_m) out << "description = " << *param.description_m << ", ";
        out << "order = " << param.order_m;

        return out;
    }
}

BOOST_FUSION_ADAPT_STRUCT(
    SDDS::parameter,
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

    template <typename Iterator>
    struct parameter_parser: qi::grammar<Iterator, parameter(), skipper<Iterator> >
    {
        parameter_parser(error_handler<Iterator> & _error_handler);

        qi::rule<Iterator, std::string(), skipper<Iterator> > string, quoted_string, units;
        qi::rule<Iterator, std::string(), skipper<Iterator> > parameter_name, parameter_units,
                parameter_description, parameter_symbol, parameter_format;
        qi::rule<Iterator, ast::datatype(), skipper<Iterator> > parameter_type;
        qi::rule<Iterator, long(), skipper<Iterator> > parameter_fixed;
        qi::rule<Iterator, parameter(), skipper<Iterator> > start;
        qi::rule<Iterator, ast::nil(), skipper<Iterator> > parameter_unsupported_pre,
                parameter_unsupported_post;
        qi::symbols<char, ast::datatype> datatype;
    };
}}
#endif /* PARAMETER_HPP_ */