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
#include "Util/SDDSParser/value_parser.hpp"

#include <array>
#include <cstddef>
#include <iostream>
#include <list>
#include <optional>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace SDDS {
    struct parameter
    {
        enum class attributes { NAME
                              , SYMBOL
                              , UNITS
                              , DESCRIPTION
                              , FORMAT_STRING
                              , TYPE
                              , FIXED_VALUE
         };

        unsigned int order_m { 0 };
        std::optional<std::string> name_m;
        std::optional<std::string> units_m;
        std::optional<std::string> description_m;
        std::optional<ast::dataType> type_m;
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
                constexpr std::array<std::pair<attributes, std::string_view>, 3> unsupportedAttributeNames = {{
                    { attributes::SYMBOL, "symbol" },
                    { attributes::FORMAT_STRING, "format_string" },
                    { attributes::FIXED_VALUE, "fixed_value" }
                }};

                for (const auto& item : unsupportedAttributeNames) {
                    if (item.first == A) {
                        std::cerr << item.second << " not supported yet" << std::endl;
                        return false;
                    }
                }

                return true;
            }
        };

        bool parse(std::string_view input, std::size_t& pos)
        {
            if (!type_m) {
                return false;
            }

            parser::ValueParser parser(input, pos);
            switch(*this->type_m) {
            case ast::dataType::FLOAT:
            {
                                float f = 0.0f;
                if (parser.parseFloat(f)) {
                    this->value_m = f;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::dataType::DOUBLE:
            {
                double d = 0.0;
                if (parser.parseDouble(d)) {
                    this->value_m = d;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::dataType::SHORT:
            {
                short s = 0;
                if (parser.parseShort(s)) {
                    this->value_m = s;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::dataType::LONG:
            {
                long l = 0;
                if (parser.parseLong(l)) {
                    this->value_m = l;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::dataType::CHARACTER:
            {
                char c = 0;
                if (parser.parseChar(c)) {
                    this->value_m = c;
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::dataType::STRING:
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
        if (param.units_m) out << "units = " << *param.units_m << ", ";
        if (param.description_m) out << "description = " << *param.description_m << ", ";
        out << "order = " << param.order_m;

        return out;
    }
}

#endif /* PARAMETER_HPP_ */