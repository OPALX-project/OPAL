//
// Class NativeHashGenerator
//   Generates a hash name.
//   Concatenates and hashes a vector of strings.
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
#include <vector>
#include <string>
#include <sstream>
#include <functional>
#include <algorithm>

class NativeHashGenerator {

public:

    static std::string generate(std::vector<std::string> arguments,
                                size_t world_pid = 0) {

        std::string hash_input = "";

        for (const std::string &arg: arguments ) {
            hash_input += arg;
        }

        hash_input += "_" + std::to_string(world_pid);

        std::hash<std::string> hashFunction;
        size_t hash_value = hashFunction(hash_input);

        std::ostringstream hash_str;
        hash_str << std::hex << hash_value;

        reverse(hash_input.begin(), hash_input.end());
        hash_value = hashFunction(hash_input);
        hash_str << hash_value;

        return hash_str.str();
    }

private:

    NativeHashGenerator() {}
    ~NativeHashGenerator() {}

};