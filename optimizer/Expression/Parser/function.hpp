//
// Namespace function
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
#ifndef FUNCTION_HPP
#define FUNCTION_HPP

#include <functional>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

namespace client {
    namespace function {
        typedef std::variant<double, bool, std::string> argument_t;

        typedef std::vector<argument_t> arguments_t;

        typedef std::function<std::tuple<double, bool> (arguments_t)> type;

        typedef std::pair<std::string, type> named_t;
    }
}

#endif