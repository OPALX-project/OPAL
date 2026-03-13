//
// Struct include
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
#ifndef INCLUDE_HPP_
#define INCLUDE_HPP_

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/error_handler.hpp"
#include "Util/SDDSParser/skipper.hpp"

#include <list>

namespace SDDS {
    struct include
    {
        enum attributes { FILENAME
                        , INCLUDE
        };

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply()
            {
                std::string attributeString;
                switch(A)
                {
                case FILENAME:
                    attributeString = "filename";
                    break;
                case INCLUDE:
                    attributeString = "include";
                    break;
                default:
                    return true;
                }
                std::cerr << attributeString << " not supported yet" << std::endl;
                return false;
            }
        };
    };

    struct includeList: std::list<include> {};

    inline std::ostream& operator<<(std::ostream& out, const include& ) {
        return out;
    }
}

#endif /* INCLUDE_HPP_ */
