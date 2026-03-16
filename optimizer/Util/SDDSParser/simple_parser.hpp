//
// Class SimpleParser
//   Simple recursive descent parser for SDDS format.
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
#ifndef SIMPLE_PARSER_HPP_
#define SIMPLE_PARSER_HPP_

#include "Util/SDDSParser/array.hpp"
#include "Util/SDDSParser/associate.hpp"
#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/column.hpp"
#include "Util/SDDSParser/data.hpp"
#include "Util/SDDSParser/description.hpp"
#include "Util/SDDSParser/file.hpp"
#include "Util/SDDSParser/include.hpp"
#include "Util/SDDSParser/parameter.hpp"
#include "Util/SDDSParser/version.hpp"

#include <cctype>
#include <stdexcept>
#include <string>
#include <string_view>

namespace SDDS {
    namespace parser {

        class SimpleParser {
        public:
            SimpleParser(std::string_view input)
                : input_(input), pos_(0), dataStartPos_(std::string::npos) {}

            file parse() {
                file result;
                dataStartPos_ = std::string::npos;
                skipWSAndComments();
                
                // Parse version
                result.sddsVersion_m = parseVersion();
                skipWSAndComments();

                // Parse optional description
                if (match("&description", true)) {
                    result.sddsDescription_m = parseDescription();
                    skipWSAndComments();
                }

                // Parse parameters, columns, associates, arrays, includes
                while (pos_ < input_.size()) {
                    skipWSAndComments();
                    if (match("&parameter", true)) {
                        result.sddsParameters_m.push_back(parseParameter());
                    } else if (match("&column", true)) {
                        result.sddsColumns_m.push_back(parseColumn());
                    } else if (match("&associate", true)) {
                        result.sddsAssociates_m.push_back(parseAssociate());
                    } else if (match("&array", true)) {
                        result.sddsArrays_m.push_back(parseArray());
                    } else if (match("&include", true)) {
                        result.sddsIncludes_m.push_back(parseInclude());
                    } else if (match("&data", true)) {
                        result.sddsData_m = parseData();
                        break;
                    } else {
                        break;
                    }
                    skipWSAndComments();
                }

                return result;
            }

            size_t dataStartPos() const {
                return dataStartPos_;
            }

        private:
            std::string_view input_;
            size_t pos_;
            size_t dataStartPos_;

            void skipWSAndComments() {
                while (pos_ < input_.size()) {
                    if (std::isspace(static_cast<unsigned char>(input_[pos_]))) {
                        pos_++;
                    } else if (input_[pos_] == '!') {
                        pos_++;
                        while (pos_ < input_.size() && input_[pos_] != '\n') pos_++;
                    } else if (input_[pos_] == ',') {
                        pos_++;
                    } else {
                        break;
                    }
                }
            }

            bool match(std::string_view str, bool ignoreCase = false) {
                skipWSAndComments();

                if (pos_ + str.size() > input_.size())
                    return false;

                const char* data = input_.data() + pos_;
                if (!ignoreCase) {
                    if (std::char_traits<char>::compare(data, str.data(), str.size()) == 0) {
                        pos_ += str.size();
                        return true;
                    }
                } else {
                    for (size_t i = 0; i < str.size(); ++i) {
                        if (std::tolower(static_cast<unsigned char>(data[i])) !=
                            std::tolower(static_cast<unsigned char>(str[i]))) {
                            return false;
                        }
                    }
                    pos_ += str.size();
                    return true;
                }

                return false;
            }

            bool match(char c) {
                skipWSAndComments();
                if (pos_ < input_.size() && input_[pos_] == c) {
                    pos_++;
                    return true;
                }
                return false;
            }

            void expect(std::string_view str) {
                if (!match(str, true)) {  // case-insensitive expect
                    throw std::runtime_error("Expected: " + std::string(str));
                }
            }

            void expect(char c) {
                if (!match(c)) {
                    throw std::runtime_error(std::string("Expected: ") + c);
                }
            }

            std::string parseQuotedString() {
                if (pos_ >= input_.size() || input_[pos_] != '"') {
                    return "";
                }
                pos_++; // skip opening quote
                std::string result;
                while (pos_ < input_.size() && input_[pos_] != '"') {
                    if (input_[pos_] == '\\' && pos_ + 1 < input_.size()) {
                        pos_++;
                        result += input_[pos_];
                    } else {
                        result += input_[pos_];
                    }
                    pos_++;
                }
                if (pos_ < input_.size() && input_[pos_] == '"') {
                    pos_++; // skip closing quote
                }
                return result;
            }

            std::string parseIdentifier() {
                skipWSAndComments();
                size_t start = pos_;
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
                return std::string(input_.substr(start, pos_ - start));
            }

            std::string parseString() {
                skipWSAndComments();
                if (pos_ < input_.size() && input_[pos_] == '"') {
                    return parseQuotedString();
                }
                return parseIdentifier();
            }

            version parseVersion() {
                version v;
                expect("SDDS");
                skipWSAndComments();
                std::string numStr;
                while (pos_ < input_.size() && std::isdigit(static_cast<unsigned char>(input_[pos_]))) {
                    numStr += input_[pos_];
                    pos_++;
                }
                if (numStr.empty()) {
                    throw std::runtime_error("Expected version number");
                }
                v.layoutVersion_m = static_cast<short>(std::stoi(numStr));
                return v;
            }

            description parseDescription() {
                description desc;

                while (true) {
                    skipWSAndComments();
                    if (match("&end", true)) {
                        break;
                    } else if (match("name", true)) {
                        expect('='); desc.text_m = parseString();
                    } else if (match("text", true)) {
                        expect('='); desc.text_m = parseString();
                    } else if (match("contents", true)) {
                        expect('='); desc.content_m = parseString();
                    } else {
                        while (pos_ < input_.size() && input_[pos_] != ',' && input_[pos_] != '\n') {
                            pos_++;
                        }
                        if (match(',')) continue;
                    }
                }

                return desc;
            }

            parameter parseParameter() {
                parameter param;

                while (pos_ < input_.size()) {
                    skipWSAndComments();

                    if (match("&end", true)) {
                        break;
                    } else if (match("name", true)) {
                        expect('='); param.name_m = parseString();
                    } else if (match("type", true)) {
                        expect('=');
                        std::string typeStr = parseString();
                        auto type = ast::parseDataType(typeStr);
                        if (!type) {
                            throw std::runtime_error("Unknown parameter type: " + typeStr);
                        }
                        param.type_m = *type;
                    } else if (match("units", true)) {
                        expect('='); param.units_m = parseString();
                    } else if (match("description", true)) {
                        expect('='); param.description_m = parseString();
                    } else {
                        while (pos_ < input_.size() && input_[pos_] != ',' && input_[pos_] != '\n') {
                            pos_++;
                        }
                        if (match(',')) continue;
                    }

                    skipWSAndComments();
                }

                if (!param.checkMandatories()) {
                    throw std::runtime_error("Parameter declaration missing mandatory attributes");
                }

                return param;
            }

            column parseColumn() {
                column col;

                while (true) {
                    skipWSAndComments();
                    if (match("&end", true)) break;
                    else if (match("name", true)) { expect('='); col.name_m = parseString(); }
                    else if (match("type", true)) {
                        expect('=');
                        std::string typeStr = parseString();
                        auto type = ast::parseDataType(typeStr);
                        if (!type) {
                            throw std::runtime_error("Unknown column type: " + typeStr);
                        }
                        col.type_m = *type;
                    }
                    else if (match("units", true)) { expect('='); col.units_m = parseString(); }
                    else if (match("description", true)) { expect('='); col.description_m = parseString(); }
                    else { while (pos_ < input_.size() && input_[pos_] != ',' && input_[pos_] != '\n') pos_++; if (match(',')) continue; }
                }

                if (!col.checkMandatories()) {
                    throw std::runtime_error("Column declaration missing mandatory attributes");
                }

                return col;
            }

            associate parseAssociate() {
                associate assoc;

                while (true) {
                    skipWSAndComments();
                    if (match("&end", true)) break;
                    parseString();
                    skipWSAndComments();
                    if (match(',')) continue;
                }

                return assoc;
            }

            array parseArray() {
                array arr;

                while (true) {
                    skipWSAndComments();
                    if (match("&end", true)) break;
                    parseString();
                    skipWSAndComments();
                    if (match(',')) continue;
                }

                return arr;
            }

            include parseInclude() {
                include inc;

                while (true) {
                    skipWSAndComments();
                    if (match("&end", true)) break;
                    parseString();
                    skipWSAndComments();
                    if (match(',')) continue;
                }

                return inc;
            }

            data parseData() {
                data d;
                skipWSAndComments();

                bool modeSet = false;
                while (pos_ < input_.size()) {
                    skipWSAndComments();
                    if (match("&end", true)) {
                        break;
                    }

                    auto key = parseIdentifier();
                    if (key.empty()) {
                        throw std::runtime_error("Malformed &data entry");
                    }
                    expect('=');

                    if (ast::equalsIgnoreCase(key, "mode")) {
                        auto modeStr = parseString();
                        if (ast::equalsIgnoreCase(modeStr, "ascii")) {
                            d.mode_m = ast::dataMode::ASCII;
                            modeSet = true;
                        } else if (ast::equalsIgnoreCase(modeStr, "binary")) {
                            d.mode_m = ast::dataMode::BINARY;
                            modeSet = true;
                        } else {
                            throw std::runtime_error("Unknown data mode: " + modeStr);
                        }
                    } else {
                        // Parse and ignore unsupported attributes while preserving cursor position.
                        parseString();
                    }

                    skipWSAndComments();
                    if (match(',')) {
                        continue;
                    }
                }

                if (!modeSet) {
                    throw std::runtime_error("Missing mode attribute in &data");
                }
                skipWSAndComments();
                dataStartPos_ = pos_;

                return d;
            }

        };

    } // namespace parser
} // namespace SDDS

#endif /* SIMPLE_PARSER_HPP_ */
