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

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/error_handler.hpp"
#include "Util/SDDSParser/skipper.hpp"
#include "Util/SDDSParser/value_parser.hpp"

#include <iostream>
#include <optional>
#include <string>
#include <vector>

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
        std::optional<std::string> name_m;
        std::optional<std::string> units_m;
        std::optional<std::string> description_m;
        std::optional<ast::datatype> type_m;
        ast::columnData_t values_m;
        static unsigned int count_m;

        bool checkMandatories() const {
            return name_m && type_m;
        }

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply() {
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

        bool parse(const std::string& input, size_t& pos) {
            parser::ValueParser parser(input, pos);
            switch(*this->type_m) {
            case ast::FLOAT:
            {
                float f = 0.0f;
                if (parser.parseFloat(f)) {
                    this->values_m.push_back(f);
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::DOUBLE:
            {
                double d = 0.0;
                if (parser.parseDouble(d)) {
                    this->values_m.push_back(d);
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::SHORT:
            {
                short s = 0;
                if (parser.parseShort(s)) {
                    this->values_m.push_back(s);
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::LONG:
            {
                long l = 0;
                if (parser.parseLong(l)) {
                    this->values_m.push_back(l);
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::CHARACTER:
            {
                char c = '\0';
                if (parser.parseChar(c)) {
                    this->values_m.push_back(c);
                    pos = parser.getPosition();
                    return true;
                }
                break;
            }
            case ast::STRING:
            {
                std::string s("");
                if (parser.parseString(s)) {
                    this->values_m.push_back(s);
                    pos = parser.getPosition();
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

#endif /* COLUMN_HPP_ */