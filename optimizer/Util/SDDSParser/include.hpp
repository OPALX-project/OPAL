//
// Struct include
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
#ifndef INCLUDE_HPP_
#define INCLUDE_HPP_

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/error_handler.hpp"

#include <array>
#include <iostream>
#include <list>
#include <string>
#include <string_view>

namespace SDDS {
    struct include
    {
        enum class attributes { FILENAME
                              , INCLUDE
        };

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply()
            {
                constexpr std::array<std::pair<attributes, std::string_view>, 2> unsupportedAttributeNames = {{
                    { attributes::FILENAME, "filename" },
                    { attributes::INCLUDE, "include" }
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

    struct includeList: std::list<include> {};

    inline std::ostream& operator<<(std::ostream& out, const include& ) {
        return out;
    }
}

#endif /* INCLUDE_HPP_ */
