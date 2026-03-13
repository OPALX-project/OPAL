//
// Struct description
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
#ifndef DESCRIPTION_HPP_
#define DESCRIPTION_HPP_

#include "Util/SDDSParser/ast.hpp"
#include "Util/SDDSParser/error_handler.hpp"
#include "Util/SDDSParser/skipper.hpp"

#include <optional>
#include <ostream>
#include <string>

namespace SDDS {
    struct description
    {
        std::optional<std::string> text_m;
        std::optional<std::string> content_m;
    };

    inline std::ostream& operator<<(std::ostream& out, const description& desc) {
        if (desc.text_m) out << "text = " << *desc.text_m << ", ";
        if (desc.content_m) out << "content = " << *desc.content_m;
        return out;
    }
}

#endif /* DESCRIPTION_HPP_ */
