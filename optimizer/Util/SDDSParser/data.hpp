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
#include "Util/SDDSParser/error_handler.hpp"
#include "Util/SDDSParser/skipper.hpp"

#include <iostream>
#include <string>

namespace SDDS {
    struct data
    {
        enum attributes { MODE
                        , LINES_PER_ROW
                        , NO_ROW_COUNT
                        , FIXED_ROW_COUNT
                        , ADDITIONAL_HEADER_LINES
                        , COLUMN_MAJOR_ORDER
                        , ENDIAN
        };

        ast::datamode mode_m;
        long numberRows_m;

        bool isASCII() const {
            if (mode_m == ast::BINARY) {
                std::cerr << "can't handle binary data yet" << std::endl;
                return false;
            }
            return true;
        }

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply() {
                std::string attributeString;
                switch(A)
                {
                case LINES_PER_ROW:
                    attributeString = "lines_per_row";
                    break;
                case NO_ROW_COUNT:
                    attributeString = "no_row_count";
                    break;
                case FIXED_ROW_COUNT:
                    attributeString = "fixed_row_count";
                    break;
                case ADDITIONAL_HEADER_LINES:
                    attributeString = "additional_header_lines";
                    break;
                case COLUMN_MAJOR_ORDER:
                    attributeString = "column_major_order";
                    break;
                case ENDIAN:
                    attributeString = "endian";
                    break;
                default:
                    return true;
                }
                std::cerr << attributeString << " not supported yet" << std::endl;
                return false;
            }
        };
    };

    inline std::ostream& operator<<(std::ostream& out, const data& data_) {
        out << "mode = " << data_.mode_m;
        return out;
    }
}

// Data parsing is now handled by SimpleParser
#endif /* DATA_HPP_ */
