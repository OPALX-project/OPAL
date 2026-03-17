//
// Class SDDSParser
//   This class writes column entries of SDDS files.
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
#include "Util/SDDSParser.h"
#include "Util/SDDSParser/simple_parser.hpp"
#include "Util/SDDSParser/skipper.hpp"

#include <cctype>
#include <cstddef>
#include <exception>
#include <fstream>
#include <string>
#include <string_view>

SDDS::SDDSParser::SDDSParser():
    sddsFileName_m("")
{ }

SDDS::SDDSParser::SDDSParser(const std::string& input):
    sddsFileName_m(input)
{ }

void SDDS::SDDSParser::setInput(const std::string& input) {
    sddsFileName_m = input;
}

SDDS::file SDDS::SDDSParser::run() {
    sddsData_m.clear();
    paramNameToID_m.clear();
    columnNameToID_m.clear();

    std::string contents = readFile();

    // Use simple parser instead of boost::spirit
    parser::SimpleParser parser(contents);
    try {
        sddsData_m = parser.parse();
    } catch (const std::exception& e) {
        throw SDDSParserException("SDDSParser::run",
                                  std::string("could not parse SDDS file: ") + e.what());
    }

    // Parse parameter values from the data section, not from file start.
    std::size_t pos = parser.dataStartPos();
    const std::size_t len = contents.size();
    if (pos == std::string::npos) {
        pos = contents.length();
    }

    auto skipWSAndComments = [&]() {
        parser::skipper::skipWhitespaceAndComments(contents, pos);
    };

    for (auto& param : sddsData_m.sddsParameters_m) {
        if (!param.parse(contents, pos)) {
            throw SDDSParserException("SDDSParser::run",
                                      "could not parse parameter value");
        }
        skipWSAndComments();
    }

    // Parse column values row-by-row. A row is only accepted if all columns parse.
    while (pos < len) {
        std::size_t rowStart = pos;
        std::size_t parsedColumns = 0;

        for (auto& col : sddsData_m.sddsColumns_m) {
            if (!col.parse(contents, pos)) {
                break;
            }
            parsedColumns++;

            skipWSAndComments();
        }

        if (parsedColumns == 0) {
            pos = rowStart;
            break;
        }

        if (parsedColumns != sddsData_m.sddsColumns_m.size()) {
            throw SDDSParserException("SDDSParser::run",
                                      "could not parse a complete data row");
        }
    }

    paramNameToID_m.reserve(sddsData_m.sddsParameters_m.size());
    std::size_t param_order = 0;
    for (const SDDS::parameter &param: sddsData_m.sddsParameters_m) {
        const std::string name = normalizeKey(*param.name_m);
        paramNameToID_m.emplace(name, param_order);
        ++ param_order;
    }

    columnNameToID_m.reserve(sddsData_m.sddsColumns_m.size());
    std::size_t col_order = 0;
    for (const SDDS::column &col: sddsData_m.sddsColumns_m) {
        const std::string name = normalizeKey(*col.name_m);
        columnNameToID_m.emplace(name, col_order);
        ++ col_order;
    }

    return sddsData_m;
}

std::string SDDS::SDDSParser::readFile() {
    std::ifstream in(sddsFileName_m, std::ios::binary);

    if (in) {
        in.seekg(0, std::ios::end);
        const std::streamsize fileSize = in.tellg();
        in.seekg(0, std::ios::beg);

        if (fileSize < 0) {
            throw SDDSParserException("SDDSParser::readSDDSFile",
                                      "could not determine size of file '" + sddsFileName_m + "'");
        }

        std::string contents(static_cast<std::size_t>(fileSize), '\0');
        if (!contents.empty()) {
            in.read(contents.data(), fileSize);
        }

        return contents;
    }

    throw SDDSParserException("SDDSParser::readSDDSFile",
                              "could not open file '" + sddsFileName_m + "'");

    return std::string("");
}

const SDDS::ast::columnData_t&
SDDS::SDDSParser::getColumnData(const std::string& columnName) const {
    const int idx = getColumnIndex(columnName);
    return sddsData_m.sddsColumns_m[idx].values_m;
}

int SDDS::SDDSParser::getColumnIndex(const std::string& col_name) const {
    const std::string normalizedName = normalizeKey(col_name);
    const auto it = columnNameToID_m.find(normalizedName);
    if (it != columnNameToID_m.end()) {
        return static_cast<int>(it->second);
    }
    throw SDDSParserException("SDDSParser::getColumnIndex",
                              "could not find column '" + col_name + "'");
}

std::string SDDS::SDDSParser::normalizeKey(std::string_view value) {
    std::size_t firstUppercase = std::string_view::npos;
    for (std::size_t index = 0; index < value.size(); ++index) {
        if (std::isupper(static_cast<unsigned char>(value[index]))) {
            firstUppercase = index;
            break;
        }
    }

    if (firstUppercase == std::string_view::npos) {
        return std::string(value);
    }

    std::string lower(value);
    for (std::size_t index = firstUppercase; index < lower.size(); ++index) {
        lower[index] = static_cast<char>(std::tolower(static_cast<unsigned char>(lower[index])));
    }

    return lower;
}
