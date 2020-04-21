//
// Class CmdArguments
//   Parsing command line arguments
//
//   In order to have a flexible framework, each component implementation gets
//   access to all command line arguments.
//   All command line options have the form:
//       --name=value
//   Spaces before and after the "=" will be trimmed.
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
#include "Util/CmdArguments.h"

void CmdArguments::addArguments(int argc, char **argv) {

    for(int i=1; i<argc; i++) {
        std::string arg = argv[i];
        std::string name, value;
        this->split(name, value, arg);
        arguments_.insert(std::pair<std::string, std::string>(name, value));
    }
}

void CmdArguments::split(std::string &name,
                         std::string &value, std::string arg) {

    size_t pos = arg.find("=");
    //strip leading '--' and '='
    name = arg.substr(2, pos - 2);
    value = arg.substr(pos + 1);

    boost::trim(name);
    boost::trim(value);
}

char** CmdArguments::getArguments() const {
    const unsigned int size = arguments_.size();
    char** args = new char*[2 * size];

    unsigned int i = 0;
    auto it = arguments_.cbegin();
    const auto end = arguments_.cend();
    for (; it != end; ++ it) {
        const std::string &key = it->first;
        char* argname = new char[key.length() + 1];
        strcpy(argname, key.c_str());
        args[i ++] = argname;

        const std::string &value = it->second;
        char* argvalue = new char[value.length() + 1];
        strcpy(argvalue, value.c_str());
        args[i ++] = argvalue;
    }

    return args;
}