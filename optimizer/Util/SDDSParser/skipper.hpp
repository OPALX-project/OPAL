//
// Namespace skipper
//
// Copyright (c) 2026, Paul Scherrer Institute, Villigen PSI, Switzerland
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
#ifndef SKIPPER_HPP_
#define SKIPPER_HPP_

#include <cctype>
#include <cstddef>
#include <string_view>

namespace SDDS {
    namespace parser {
        namespace skipper {

            inline void skipWhitespaceAndComments(std::string_view input, std::size_t& pos) {
                while (pos < input.size()) {
                    if (std::isspace(static_cast<unsigned char>(input[pos]))) {
                        pos++;
                    } else if (input[pos] == '!') {
                        // Skip comment to end of line.
                        pos++;
                        while (pos < input.size() && input[pos] != '\n') {
                            pos++;
                        }
                    } else if (input[pos] == ',') {
                        // Keep comma behavior aligned with legacy parser semantics.
                        pos++;
                    } else {
                        break;
                    }
                }
            }

        } // namespace utils
    } // namespace parser
} // namespace SDDS

#endif /* SKIPPER_HPP_ */