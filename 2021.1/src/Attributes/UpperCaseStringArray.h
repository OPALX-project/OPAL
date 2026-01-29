//
// Class UpperCaseStringArray
//   This class is used to parse attributes of type array of strings that should
//   be all upper case, i.e. it returns an array of upper case strings.
//
// Copyright (c) 2020, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_UpperCaseStringArray_HH
#define OPAL_UpperCaseStringArray_HH

#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/AttributeHandler.h"
#include "AbstractObjects/Expressions.h"


namespace Attributes {

    /// Parser for an attribute of type string array.
    class UpperCaseStringArray: public AttributeHandler {

    public:

        /// Constructor.
        //  Assign attribute name and help string.
        UpperCaseStringArray(const std::string &name, const std::string &help);

        virtual ~UpperCaseStringArray();

        /// Return attribute type string ``string array''.
        virtual const std::string &getType() const;

        /// Parse the attribute.
        virtual void parse(Attribute &, Statement &, bool) const;

        /// Parse a component of the array.
        //  Identified by its index.
        virtual void parseComponent(Attribute &, Statement &, bool, int) const;

    private:

        // Not implemented.
        UpperCaseStringArray();
        UpperCaseStringArray(const UpperCaseStringArray &);
        void operator=(const UpperCaseStringArray &);
    };

};

#endif
