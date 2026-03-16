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
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
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

        constexpr std::array<std::pair<std::string_view, dataType>, 6> dataTypeTable = {{
            { "float", dataType::FLOAT },
            { "double", dataType::DOUBLE },
            { "short", dataType::SHORT },
            { "long", dataType::LONG },
            { "character", dataType::CHARACTER },
            { "string", dataType::STRING }
        }};

        constexpr char toLowerAscii(char c) {
            return (c >= 'A' && c <= 'Z') ? static_cast<char>(c + ('a' - 'A')) : c;
        }

        constexpr bool equalsIgnoreCase(std::string_view lhs, std::string_view rhs) {
            if (lhs.size() != rhs.size()) {
                return false;
            }
            for (std::size_t i = 0; i < lhs.size(); ++i) {
                if (toLowerAscii(lhs[i]) != toLowerAscii(rhs[i])) {
                    return false;
                }
            }
            return true;
        }

        constexpr std::optional<dataType> parseDataType(std::string_view typeName) {
            for (const auto& [name, type] : dataTypeTable) {
                if (equalsIgnoreCase(typeName, name)) {
                    return type;
                }
            }
            return std::nullopt;
        }

        constexpr std::string_view getDataTypeString(dataType type) {
            for (const auto& [name, value] : dataTypeTable) {
                if (value == type) {
                    return name;
                }
            }

            return "unknown";
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