//
// Class StringConstant
//   The STRING CONSTANT definition.
//
// Copyright (c) 2000 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_StringConstant_HH
#define OPAL_StringConstant_HH

#include "AbstractObjects/ValueDefinition.h"


class StringConstant: public ValueDefinition {

public:

    /// Exemplar constructor.
    StringConstant();

    virtual ~StringConstant();

    /// Test if object can be replaced.
    //  True, if [b]rhs[/b] is a string constant.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual StringConstant *clone(const std::string &name);

    /// Print the constant.
    virtual void print(std::ostream &) const;

    /// Print its value
    virtual void printValue(std::ostream &os) const;

    /// Return value.
    virtual std::string getString() const;

private:

    // Not implemented.
    StringConstant(const StringConstant &);
    void operator=(const StringConstant &);

    // Clone constructor.
    StringConstant(const std::string &name, StringConstant *parent);
    StringConstant(const std::string &name, StringConstant *parent, const std::string &value);
};

#endif // OPAL_StringConstant_HH