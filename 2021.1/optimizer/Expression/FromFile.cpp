//
// Struct FromFile
//   Simple functor that reads vector data from a file. If the file contains
//   more than one value the sum is returned.
//   \f[
//     result = \sum_{i=0}^n value_i
//   \f]
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
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
#include "Expression/FromFile.h"

/// reads a simple list of double values
void FromFile::readValues() {

    values_.clear();

    std::ifstream file;
    file.open(filename_.c_str(), std::ios::in);
    if(!file) {
        throw OptPilotException("FromFile::readValues()",
                "Error opening file " + filename_);
    }

    std::copy(std::istream_iterator<double>(file),
              std::istream_iterator<double>(),
              std::back_inserter(values_));

    file.close();
}

