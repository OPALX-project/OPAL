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
        enum class dataType { FLOAT
                            , DOUBLE
                            , SHORT
                            , LONG
                            , CHARACTER
                            , STRING };

        enum class dataMode { ASCII
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

        constexpr std::string_view getDataTypeString(dataType type) {
            constexpr std::array<std::string_view, 6> dataTypeNames = {
                "float",
                "double",
                "short",
                "long",
                "char",
                "string"
            };

            const auto index = static_cast<std::underlying_type_t<dataType>>(type);
            if (index < 0 || static_cast<std::size_t>(index) >= dataTypeNames.size()) {
                return "unknown";
            }

            return dataTypeNames[static_cast<std::size_t>(index)];
        }

        constexpr std::string_view getDataModeString(dataMode mode) {
            constexpr std::array<std::string_view, 2> dataModeNames = {
                "ascii",
                "binary"
            };

            const auto index = static_cast<std::underlying_type_t<dataMode>>(mode);
            if (index < 0 || static_cast<std::size_t>(index) >= dataModeNames.size()) {
                return "unknown";
            }

            return dataModeNames[static_cast<std::size_t>(index)];
        }

    }

    namespace parser {}
}

#endif // AST_HPP_