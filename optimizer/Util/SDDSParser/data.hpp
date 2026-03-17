//
// Struct data
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
#ifndef DATA_HPP_
#define DATA_HPP_

#include "Util/SDDSParser/ast.hpp"
#include <array>
#include <iostream>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>

namespace SDDS {
    struct data
    {
        enum class attributes { MODE
                              , LINES_PER_ROW
                              , NO_ROW_COUNT
                              , FIXED_ROW_COUNT
                              , ADDITIONAL_HEADER_LINES
                              , COLUMN_MAJOR_ORDER
                              , ENDIAN
        };

        ast::dataMode mode_m { ast::dataMode::ASCII };
        long numberRows_m { 0 };

        bool isASCII() const {
            if (mode_m == ast::dataMode::BINARY) {
                std::cerr << "can't handle binary data yet" << std::endl;
                return false;
            }
            return true;
        }

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply() {
                constexpr std::array<std::pair<attributes, std::string_view>, 6> unsupportedAttributeNames = {{
                    { attributes::LINES_PER_ROW, "lines_per_row" },
                    { attributes::NO_ROW_COUNT, "no_row_count" },
                    { attributes::FIXED_ROW_COUNT, "fixed_row_count" },
                    { attributes::ADDITIONAL_HEADER_LINES, "additional_header_lines" },
                    { attributes::COLUMN_MAJOR_ORDER, "column_major_order" },
                    { attributes::ENDIAN, "endian" }
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

    inline std::ostream& operator<<(std::ostream& out, const data& data_) {
        out << "mode = " << ast::getDataModeString(data_.mode_m);
        return out;
    }
}

// Data parsing is now handled by SimpleParser
#endif /* DATA_HPP_ */
