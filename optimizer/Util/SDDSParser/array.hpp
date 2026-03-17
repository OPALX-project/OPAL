//
// Struct array
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
#ifndef ARRAY_HPP_
#define ARRAY_HPP_

#include <array>
#include <iostream>
#include <list>
#include <ostream>
#include <string_view>
#include <utility>

namespace SDDS {
    struct array {
        enum class attributes {
            NAME,
            SYMBOL,
            UNITS,
            DESCRIPTION,
            FORMAT_STRING,
            GROUP_NAME,
            TYPE,
            FIELD_LENGTH,
            DIMENSIONS,
            ARRAY
         };

        template <attributes A>
        struct complainUnsupported {
            static bool apply() {
                constexpr std::array<std::pair<attributes, std::string_view>, 10> unsupportedAttributeNames = {{
                    { attributes::NAME, "name" },
                    { attributes::SYMBOL, "symbol" },
                    { attributes::UNITS, "units" },
                    { attributes::DESCRIPTION, "description" },
                    { attributes::FORMAT_STRING, "format_string" },
                    { attributes::GROUP_NAME, "group_name" },
                    { attributes::TYPE, "type" },
                    { attributes::FIELD_LENGTH, "field_length" },
                    { attributes::DIMENSIONS, "dimensions" },
                    { attributes::ARRAY, "array" }
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
    };

    struct arrayList : std::list<array> {};

    inline std::ostream& operator<<(std::ostream& out, const array&) {
        return out;
    }
}

#endif /* ARRAY_HPP_ */