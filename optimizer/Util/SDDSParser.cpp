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

#include <algorithm>
#include <cctype>
#include <limits>
#include <string>
#include <utility>


SDDS::SDDSParser::SDDSParser():
    sddsFileName_m("")
{ }

SDDS::SDDSParser::SDDSParser(const std::string &input):
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
    size_t pos = parser.dataStartPos();
    if (pos == std::string::npos) {
        pos = contents.length();
    }
    for (auto& param : sddsData_m.sddsParameters_m) {
        if (!param.parse(contents, pos)) {
            throw SDDSParserException("SDDSParser::run",
                                      "could not parse parameter value");
        }
        // Skip comma or whitespace
        while (pos < contents.length() && (std::isspace(static_cast<unsigned char>(contents[pos])) || contents[pos] == ',')) {
            pos++;
        }
    }

    // Parse column values row-by-row. A row is only accepted if all columns parse.
    while (pos < contents.length()) {
        size_t rowStart = pos;
        size_t parsedColumns = 0;

        for (auto& col : sddsData_m.sddsColumns_m) {
            if (!col.parse(contents, pos)) {
                break;
            }
            parsedColumns++;

            while (pos < contents.length() &&
                   (std::isspace(static_cast<unsigned char>(contents[pos])) || contents[pos] == ',')) {
                pos++;
            }
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

    unsigned int param_order = 0;
    for (const SDDS::parameter &param: sddsData_m.sddsParameters_m) {
        std::string name = *param.name_m;
        fixCaseSensitivity(name);
        paramNameToID_m.insert(std::make_pair(name,
                                              param_order));
        ++ param_order;
    }

    unsigned int col_order = 0;
    for (const SDDS::column &col: sddsData_m.sddsColumns_m) {
        std::string name = *col.name_m;
        fixCaseSensitivity(name);
        columnNameToID_m.insert(std::make_pair(name,
                                               col_order));
        ++ col_order;
    }

    return sddsData_m;
}

std::string SDDS::SDDSParser::readFile() {
    std::ifstream in(sddsFileName_m.c_str());

    if (in) {
        std::string contents;
        in.seekg(0, std::ios::end);
        contents.resize(in.tellg());
        in.seekg(0, std::ios::beg);

        in.read(&contents[0], contents.size());

        in.close();

        return contents;
    }

    throw SDDSParserException("SDDSParser::readSDDSFile",
                              "could not open file '" + sddsFileName_m + "'");

    return std::string("");
}

SDDS::ast::columnData_t SDDS::SDDSParser::getColumnData(const std::string& columnName) {
    int idx = getColumnIndex(columnName);

    return sddsData_m.sddsColumns_m[idx].values_m;
}

int SDDS::SDDSParser::getColumnIndex(std::string col_name) const {
    fixCaseSensitivity(col_name);
    auto it = columnNameToID_m.find(col_name);
    if (it != columnNameToID_m.end()) {
        return it->second;
    }

    throw SDDSParserException("SDDSParser::getColumnIndex",
                              "could not find column '" + col_name + "'");

}

//XXX use either all upper, or all lower case chars
void SDDS::SDDSParser::fixCaseSensitivity(std::string& for_string) {
    std::transform(for_string.begin(),
                   for_string.end(),
                   for_string.begin(),
                   [](unsigned char c) {
                        return std::tolower(c);
                    }
    );
}