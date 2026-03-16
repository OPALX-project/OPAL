//
// Namespace ast
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
#ifndef AST_HPP_
#define AST_HPP_

#include <array>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

namespace SDDS {
    namespace ast {
        enum class datatype { FLOAT
                            , DOUBLE
                            , SHORT
                            , LONG
                            , CHARACTER
                            , STRING };

        enum class datamode { ASCII
                            , BINARY};

        enum class endianess { BIGENDIAN
                             , LITTLEENDIAN};

        struct nil {};

        using variant_t = std::variant<float,
                           double,
                           short,
                           long,
                           char,
                           std::string>;

        using columnData_t = std::vector<variant_t>;

        constexpr std::string_view getDataTypeString(datatype type) {
            constexpr std::array<std::string_view, 6> datatypeNames = {
                "float",
                "double",
                "short",
                "long",
                "char",
                "string"
            };

            const auto index = static_cast<std::underlying_type_t<datatype>>(type);
            if (index < 0 || static_cast<std::size_t>(index) >= datatypeNames.size()) {
                return "unknown";
            }

            return datatypeNames[static_cast<std::size_t>(index)];
        }

        constexpr std::string_view getDataModeString(datamode mode) {
            constexpr std::array<std::string_view, 2> datamodeNames = {
                "ascii",
                "binary"
            };

            const auto index = static_cast<std::underlying_type_t<datamode>>(mode);
            if (index < 0 || static_cast<std::size_t>(index) >= datamodeNames.size()) {
                return "unknown";
            }

            return datamodeNames[static_cast<std::size_t>(index)];
        }

    }

    namespace parser {}
}

#endif // AST_HPP_