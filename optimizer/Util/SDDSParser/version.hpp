//
// Struct version
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
#ifndef VERSION_HPP_
#define VERSION_HPP_

#include <ostream>

namespace SDDS {
    struct version {
        short layoutVersion_m { 0 };
    };

    inline std::ostream& operator<<(std::ostream& out, const version& head) {
        out << "layout version is " << head.layoutVersion_m;
        return out;
    }
}

#endif /* VERSION_HPP_ */
