//
// Class UpperCaseString
//   This class is used to parse attributes of type string that should
//   be upper case, i.e. it returns upper case strings.
//
// Copyright (c) 2020, Christof Metzger-Kraus, Open Sourcerer
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
#ifndef OPAL_UpperCaseString_HH
#define OPAL_UpperCaseString_HH

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"


namespace Attributes {

    /// Parser for an attribute of type string.
    class UpperCaseString: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        UpperCaseString(const std::string &name, const std::string &help);

        virtual ~UpperCaseString();

        /// Return attribute type string ``string''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

    private:

        // Not implemented.
        UpperCaseString();
        UpperCaseString(const UpperCaseString &);
        void operator=(const UpperCaseString &);
    };

};

#endif