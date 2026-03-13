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

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/error_handler.hpp"
#include "Util/SDDSParser/skipper.hpp"
#include "Util/SDDSParser/value_parser.hpp"

#include <list>
#include <optional>
#include <ostream>
#include <string>

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
        std::optional<std::string> name_m;
        std::optional<std::string> units_m;
        std::optional<std::string> description_m;
        std::optional<ast::datatype> type_m;
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

        bool parse(const std::string& input, size_t& pos)
        {
            parser::ValueParser parser(input, pos);
            switch(*this->type_m) {
            case ast::FLOAT:
            {
                                float f = 0.0f;
                if (parser.parseFloat(f)) {
                    this->value_m = f;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::DOUBLE:
            {
                double d = 0.0;
                if (parser.parseDouble(d)) {
                    this->value_m = d;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::SHORT:
            {
                short s = 0;
                if (parser.parseShort(s)) {
                    this->value_m = s;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::LONG:
            {
                long l = 0;
                if (parser.parseLong(l)) {
                    this->value_m = l;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::CHARACTER:
            {
                char c = 0;
                if (parser.parseChar(c)) {
                    this->value_m = c;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::STRING:
            {
                std::string s;
                if (parser.parseStringToEol(s)) {
                    this->value_m = s;
                    pos = parser.getPosition();
                    return true;
                }
                break;
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

#endif /* PARAMETER_HPP_ */