//
// Struct file
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
#ifndef FILE_HPP_
#define FILE_HPP_

#include "Util/SDDSParser/array.hpp"
#include "Util/SDDSParser/associate.hpp"
#include "Util/SDDSParser/column.hpp"
#include "Util/SDDSParser/data.hpp"
#include "Util/SDDSParser/description.hpp"
#include "Util/SDDSParser/include.hpp"
#include "Util/SDDSParser/parameter.hpp"
#include "Util/SDDSParser/version.hpp"

#include <optional>
#include <ostream>

namespace SDDS {
    struct file {
        version sddsVersion_m;                              // 0
        std::optional<description> sddsDescription_m;       // 1
        parameterList sddsParameters_m;                     // 2
        columnList sddsColumns_m;                           // 3
        data sddsData_m;                                    // 4
        associateList sddsAssociates_m;                     // 5
        arrayList sddsArrays_m;                             // 6
        includeList sddsIncludes_m;                         // 7

        void clear() {
            sddsVersion_m = {};
            sddsDescription_m.reset();
            sddsParameters_m.clear();
            sddsColumns_m.clear();
            sddsData_m = {};
            sddsAssociates_m.clear();
            sddsArrays_m.clear();
            sddsIncludes_m.clear();
        }
    };

    inline std::ostream& operator<<(std::ostream& out, const file& fh) {
        out << fh.sddsVersion_m << "\n";
        if (fh.sddsDescription_m) out << *fh.sddsDescription_m << "\n";
        for (const parameter& param: fh.sddsParameters_m) {
            out << param << "\n";
        }
        for (const column& col: fh.sddsColumns_m) {
            out << col << "\n";
        }
        out << fh.sddsData_m << "\n";

//        if (fh.sddsAssociates_m) out << *fh.sddsAssociates_m << "\n";
//        if (fh.sddsArrays_m) out << *fh.sddsArrays_m << "\n";
//        if (fh.sddsIncludes_m) out << *fh.sddsIncludes_m << "\n";

        out << std::endl;

        return out;
    }
}

#endif /* FILE_HPP_ */