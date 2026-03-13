//
// Struct associate
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
#ifndef ASSOCIATE_HPP_
#define ASSOCIATE_HPP_

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/error_handler.hpp"
#include "Util/SDDSParser/skipper.hpp"

#include <iostream>
#include <list>
#include <string>

namespace SDDS {
    struct associate
    {
        enum attributes { NAME
                        , FILENAME
                        , PATH
                        , DESCRIPTION
                        , CONTENTS
                        , SDDS
                        , ASSOCIATE
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
                case FILENAME:
                    attributeString = "filename";
                    break;
                case PATH:
                    attributeString = "path";
                    break;
                case DESCRIPTION:
                    attributeString = "description";
                    break;
                case CONTENTS:
                    attributeString = "contents";
                    break;
                case SDDS:
                    attributeString = "sdds";
                    break;
                case ASSOCIATE:
                    attributeString = "associate";
                    break;
                default:
                    return true;
                }
                std::cerr << attributeString << " not supported yet" << std::endl;
                return false;
            }
        };
    };

    struct associateList: std::list<associate> {};

    inline std::ostream& operator<<(std::ostream& out, const associate& ) {
        return out;
    }
}

#endif /* ASSOCIATE_HPP_ */
