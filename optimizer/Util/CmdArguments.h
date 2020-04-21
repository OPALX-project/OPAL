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
#ifndef __CMD_ARGUMENTS__
#define __CMD_ARGUMENTS__

#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <set>

#include "boost/smart_ptr.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/utility/value_init.hpp"

#include "Util/OptPilotException.h"

class CmdArguments {

public:

    CmdArguments(int argc, char **argv) {
        addArguments(argc, argv);
    }

    CmdArguments()
    {}

    ~CmdArguments()
    {}

    /** Get command line argument without specifying default value.
     *  @see getArg
    */
    template<class T>
    T getArg(const std::string name, bool isFatal = false) {
        boost::value_initialized<T> value;
        return getArg<T>(name, value, isFatal);
    }

    /** Get command line argument with default value. If the parameter is not
     *  found either an exception is thrown or a warning is written to stdout
     *  and the default value is used (@see isFatal).
     *
     *  @param[in] name of the command line argument
     *  @param[in] default_value of the parameter
     *  @param[in] isFatal throw exception if parameter is expected to be
     *             present, default is false
     *  @return argument value
     */
    template<class T>
    T getArg(const std::string name, T default_value, bool isFatal = false) {
        T value = default_value;
        try {
            value = this->arg<T>(name);
        } catch(OptPilotException &e) {
            if(isFatal) {
                std::ostringstream exe;
                exe << "Fatal: argument \"";
                exe << name;
                exe << "\" not found!";
                throw OptPilotException("CmdArguments::getArg", exe.str());
            } else {
                if(warned_.count(name) == 0) {
                    std::ostringstream warn;
                    warn << "\033[01;35m";
                    warn << "Warning: argument \"";
                    warn << name;
                    warn << "\" not found! Using default value (";
                    warn << default_value;
                    warn << ").";
                    warn << "\e[0m" << std::endl;
                    std::cout << warn.str() << std::flush;
                    warned_.insert(name);
                }
            }
        }

        return value;
    }

    template <class T>
    void addArgument(const std::string &name, const T &value) {
        std::ostringstream oss;
        oss << value;

        if (arguments_.find(name) != arguments_.end()) {
            throw OptPilotException("CmdArguments::addArgument",
                                    "Argument '" + name + "' exists");
        }

        arguments_.insert(std::make_pair(name, oss.str()));
    }

    template <class T>
    void replaceArgument(const std::string &name, const T &value) {
        std::ostringstream oss;
        oss << value;

        if (arguments_.find(name) == arguments_.end()) {
            arguments_.insert(std::make_pair(name, oss.str()));
            return;
        }

        arguments_.at(name) = oss.str();
    }

    char** getArguments() const;
    unsigned int getNumArguments() const {
        return arguments_.size();
    }
    //template<>
    //size_t getArg<size_t>(const std::string name, bool isFatal = false) {
        //return getArg<size_t>(name, 0, isFatal);
    //}

    //template<>
    //int getArg<int>(const std::string name, bool isFatal = false) {
        //return getArg<int>(name, 0, isFatal);
    //}

    //template<>
    //unsigned int getArg<unsigned int>(const std::string name, bool isFatal = false) {
        //return getArg<unsigned int>(name, 0, isFatal);
    //}

    //template<>
    //std::string getArg<std::string>(const std::string name, bool isFatal = false) {
        //return getArg<std::string>(name, "", isFatal);
    //}

private:

    /// container for storing command line options
    std::map<std::string, std::string> arguments_;

    std::set<std::string> warned_;

    /// parse user arguments
    void addArguments(int argc, char **argv);

    /// helper to split string
    void split(std::string &name, std::string &value, std::string arg);

    /// tries to retrieve command line parameter.
    /// @throws OptPilotException if parameter was not found.
    template<class T>
    T arg(const std::string name);

};

typedef boost::shared_ptr<CmdArguments> CmdArguments_t;

template<class T>
inline T CmdArguments::arg(const std::string name) {
    T t;
    std::map<std::string, std::string>::iterator it = arguments_.find(name);
    if(it != arguments_.end()) {
        std::istringstream iss(arguments_[name]);
        iss >> t;
        return t;
    } else {
        throw OptPilotException("CmdArguments::getArg", "argument not found!");
    }
}

template<>
inline std::string CmdArguments::arg<std::string>(const std::string name) {
    std::map<std::string, std::string>::iterator it = arguments_.find(name);
    if(it != arguments_.end()) {
        return arguments_[name];
    } else {
        throw OptPilotException("CmdArguments::getArg", "argument not found!");
    }
}

#endif