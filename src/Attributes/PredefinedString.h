//
// Class PredefinedString
//   This class is used to parse attributes of type string that should
//   be contained in set of predefined strings.
//
// Copyright (c) 2021, Christof Metzger-Kraus, Open Sourcerer
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
#ifndef OPAL_PredefinedString_HH
#define OPAL_PredefinedString_HH

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"

#include <set>

namespace Attributes {

    /// Parser for an attribute of type string.
    class PredefinedString: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        PredefinedString(const std::string& name,
                         const std::string& help,
                         const std::initializer_list<std::string>& predefinedStrings,
                         const std::string& defaultValue = "_HAS_NO_DEFAULT_");

        virtual ~PredefinedString();

        /// Return attribute type string ``string''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        PredefinedString();
        PredefinedString(const PredefinedString &);
        void operator=(const PredefinedString &);

        std::set<std::string> predefinedStrings_m;
    };

};

#endif