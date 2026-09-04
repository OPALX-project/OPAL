//
// Class FromFile
//   This class parses a file that contains design variable values.
//   Each column belongs to a design variable.
//   The first line is considered as headerLine and consists of the
//   design variable name. The name has to agree with the string
//   in the input file.
//
// Copyright (c) 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#include "Sample/FromFile.h"

#include "Utilities/OpalException.h"

#include <algorithm>
#include <fstream>
#include <iterator>
#include <sstream>

FromFile::FromFile(const std::string& filename, const std::string& dvarName, std::size_t modulo)
    : mod_m(modulo), filename_m(filename), dvarName_m(dvarName), globalSize_m(0)
{
    // Open the input file
    std::ifstream in(filename_m);
    if (!in.is_open()) {
        throw OpalException("FromFile()",
                            "Couldn't open file '" + filename_m + "'.");
    }

    // Count number of non-empty lines
    std::string line;
    std::size_t validLines = 0;
    while (std::getline(in, line)) {
        if (!line.empty()) ++validLines;
    }
    if (validLines == 0) {
        throw OpalException("FromFile()", "Empty file '" + filename_m + "'.");
    }

    // Store the number of data rows (including header)
    globalSize_m = static_cast<unsigned int>(validLines);

    in.close();
}

void FromFile::create(std::shared_ptr<SampleIndividual>& ind, std::size_t i) {
    ind->genes[i] = getNext(ind->id);
}

void FromFile::allocate(const CmdArguments_t& /*args*/, const Comm::Bundle_t& /*comm*/) {
    std::ifstream file(filename_m);
    if (!file.is_open()) {
        throw OpalException("FromFile::allocate",
                "Couldn't open file '" + filename_m + "'.");
    }

    // Read and parse the header line to identify the variable names
    std::string headerLine;
    if (!std::getline(file, headerLine)) {
        throw OpalException("FromFile::allocate",
                "Empty file or error reading header in '" + filename_m + "'.");
    }

    std::istringstream headerStream(headerLine);
    std::vector<std::string> dvars{std::istream_iterator<std::string>{headerStream},
                                   std::istream_iterator<std::string>{}};

    // Find the column index of the requested design variable
    auto it = std::find(dvars.begin(), dvars.end(), dvarName_m);
    if (it == dvars.end()) {
        throw OpalException("FromFile::allocate",
                "Couldn't find dvar '" + dvarName_m +
                "' in file '" + filename_m + "'");
    }
    const std::size_t columnIndex = std::distance(dvars.begin(), it);

    // Reserve space for data values (excluding the header)
    chain_m.reserve(globalSize_m - 1); // -1 for the header line
    std::string line;
    std::size_t lineNumber = 2; // Start at 2 because the header is line 1
    while (std::getline(file, line)) {
        if (line.empty()) {
            ++lineNumber;
            continue;
        }

        std::istringstream lineStream(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{lineStream},
                                        std::istream_iterator<std::string>{}};

        if (columnIndex >= tokens.size()) {
            throw OpalException("FromFile::allocate",
                    "Line " + std::to_string(lineNumber) + " in file '" + filename_m +
                    "' has fewer columns (" + std::to_string(tokens.size()) +
                    ") than expected (index " + std::to_string(columnIndex) + ").");
        }

        try {
            chain_m.push_back(std::stod(tokens[columnIndex]));
        } catch (const std::exception& e) {
            throw OpalException("FromFile::allocate",
                    "Conversion error at line " + std::to_string(lineNumber) +
                    " in file '" + filename_m + "': " + e.what());
        }

        ++lineNumber;
    }

    file.close();
}

double FromFile::getNext(unsigned int id) {
    if (globalSize_m == 0 || chain_m.empty()) {
        throw OpalException("FromFile::getNext",
            "No valid numeric data found for DVAR '" + dvarName_m +
            "' in file '" + filename_m + "'.");
    }

    std::size_t idx = (id / mod_m) % chain_m.size();
    return chain_m[idx];
}

unsigned int FromFile::getSize() const {
    return globalSize_m;
}
