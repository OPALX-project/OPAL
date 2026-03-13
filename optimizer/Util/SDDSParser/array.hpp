//
// Struct array
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
#ifndef ARRAY_HPP_
#define ARRAY_HPP_

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/error_handler.hpp"
#include "Util/SDDSParser/skipper.hpp"

#include <iostream>
#include <list>
#include <string>

namespace SDDS {
    struct array
    {
        enum attributes { NAME
                         , SYMBOL
                         , UNITS
                         , DESCRIPTION
                         , FORMAT_STRING
                         , GROUP_NAME
                         , TYPE
                         , FIELD_LENGTH
                         , DIMENSIONS
                         , ARRAY
         };

        template <attributes A>
        struct complainUnsupported
        {
            static bool apply()
            {
                std::string attributeString;
                switch(A)
                {
                case NAME:
                    attributeString = "name";
                    break;
                case SYMBOL:
                    attributeString = "symbol";
                    break;
                case UNITS:
                    attributeString = "units";
                    break;
                case DESCRIPTION:
                    attributeString = "description";
                    break;
                case FORMAT_STRING:
                    attributeString = "format_string";
                    break;
                case GROUP_NAME:
                    attributeString = "group_name";
                    break;
                case TYPE:
                    attributeString = "type";
                    break;
                case FIELD_LENGTH:
                    attributeString = "field_length";
                    break;
                case DIMENSIONS:
                    attributeString = "dimensions";
                    break;
                case ARRAY:
                    attributeString = "array";
                    break;
                default:
                    return true;
                }
                std::cerr << attributeString << " not supported yet" << std::endl;
                return false;
            }
        };
    };

    struct arrayList: std::list<array> {};

    inline std::ostream& operator<<(std::ostream& out, const array& ) {
        return out;
    }
}

#endif /* ARRAY_HPP_ */
