//
// Class ValueParser
//   Simple value parser for SDDS data values.
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
#ifndef VALUE_PARSER_HPP_
#define VALUE_PARSER_HPP_

#include "Util/SDDSParser/skipper.hpp"

#include <cctype>
#include <string>
#include <string_view>
#include <type_traits>

namespace SDDS {
    namespace parser {

        class ValueParser {
        public:
            ValueParser(std::string_view input, std::size_t start_pos = 0)
                : input_(input), pos_(start_pos) {}

            void skipWSAndComments() {
                skipper::skipWhitespaceAndComments(input_, pos_);
            }

            // ----------------------- Numeric Parsers -----------------------
            template <typename T>
            bool parseNumeric(T& value) {
                skipWSAndComments();
                std::size_t start = pos_;
                bool has_dot = false;
                bool has_exp = false;

                if (pos_ < input_.size() && (input_[pos_] == '+' || input_[pos_] == '-')) {
                    pos_++;
                }

                while (pos_ < input_.size()) {
                    char c = input_[pos_];
                    if (std::isdigit(static_cast<unsigned char>(c))) {
                        pos_++;
                    } else if (c == '.' && !has_dot) {
                        has_dot = true;
                        pos_++;
                    } else if ((c == 'e' || c == 'E') && !has_exp) {
                        has_exp = true;
                        pos_++;
                        if (pos_ < input_.size() && (input_[pos_] == '+' || input_[pos_] == '-')) pos_++;
                    } else {
                        break;
                    }
                }

                if (pos_ > start) {
                    try {
                        const std::string token(input_.data() + start, pos_ - start);
                        if constexpr (std::is_same_v<T, float>) {
                            value = std::stof(token);
                        } else {
                            value = std::stod(token);
                        }
                        return true;
                    } catch (...) {
                        pos_ = start;
                        return false;
                    }
                }
                pos_ = start;
                return false;
            }

            bool parseInteger(long& value) {
                skipWSAndComments();
                std::size_t start = pos_;
                if (pos_ < input_.size() && (input_[pos_] == '+' || input_[pos_] == '-')) {
                    pos_++;
                }

                while (pos_ < input_.size() && std::isdigit(static_cast<unsigned char>(input_[pos_]))) {
                    pos_++;
                }

                if (pos_ > start) {
                    try {
                        const std::string token(input_.data() + start, pos_ - start);
                        value = std::stol(token);
                        return true;
                    } catch (...) {
                        pos_ = start;
                        return false;
                    }
                }
                pos_ = start;
                return false;
            }

            bool parseFloat(float& value) { return parseNumeric(value); }
            bool parseDouble(double& value) { return parseNumeric(value); }

            bool parseShort(short& value) {
                long tmp;
                if (!parseInteger(tmp)) return false;
                value = static_cast<short>(tmp);
                return true;
            }

            bool parseLong(long& value) { return parseInteger(value); }

            bool parseChar(char& value) {
                skipWSAndComments();
                if (pos_ < input_.size() && !std::isspace(static_cast<unsigned char>(input_[pos_]))) {
                    value = input_[pos_++];
                    return true;
                }
                return false;
            }

            // ----------------------- String Parsers -----------------------
            bool parseQuotedString(std::string& value) {
                skipWSAndComments();
                if (pos_ >= input_.size() || input_[pos_] != '"') {
                    return false;
                }
                pos_++;
                value.clear();
                while (pos_ < input_.size() && input_[pos_] != '"') {
                    if (input_[pos_] == '\\' && pos_ + 1 < input_.size()) {
                        pos_++;
                        value += input_[pos_];
                    } else {
                        value += input_[pos_];
                    }
                    pos_++;
                }
                if (pos_ < input_.size() && input_[pos_] == '"') {
                    pos_++;
                    return true;
                }
                return false;
            }

            bool parseString(std::string& value) {
                skipWSAndComments();
                if (pos_ < input_.size() && input_[pos_] == '"') {
                    return parseQuotedString(value);
                }

                std::size_t start = pos_;
                while (pos_ < input_.size()) {
                    char c = input_[pos_];
                    if (std::isalnum(static_cast<unsigned char>(c)) ||
                        c == '@' || c == '#' || c == ':' || c == '+' ||
                        c == '-' || c == '%' || c == '.' || c == '_' ||
                        c == '$' || c == '&' || c == '/') {
                        pos_++;
                    } else {
                        break;
                    }
                }
                if (pos_ > start) {
                    value.assign(input_.data() + start, pos_ - start);
                    return true;
                }
                return false;
            }

            bool parseStringToEol(std::string& value) {
                skipWSAndComments();
                if (pos_ < input_.size() && input_[pos_] == '"') {
                    return parseQuotedString(value);
                }

                std::size_t start = pos_;
                while (pos_ < input_.size() && input_[pos_] != '\n' && input_[pos_] != '\r') {
                    pos_++;
                }

                std::size_t end = pos_;
                while (end > start && (input_[end - 1] == ' ' || input_[end - 1] == '\t')) {
                    end--;
                }

                if (end > start) {
                    value.assign(input_.data() + start, end - start);
                    return true;
                }

                return false;
            }

            // ----------------------- Position Control -----------------------
            std::size_t getPosition() const { return pos_; }
            void setPosition(std::size_t pos) { pos_ = pos; }

        private:
            std::string_view input_;
            std::size_t pos_;
        };

    } // namespace parser
} // namespace SDDS

#endif /* VALUE_PARSER_HPP_ */