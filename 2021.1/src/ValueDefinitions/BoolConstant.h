//
// Class BoolConstant
//   The BOOL CONSTANT definition.
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
#ifndef OPAL_BoolConstant_HH
#define OPAL_BoolConstant_HH

#include "AbstractObjects/ValueDefinition.h"


class BoolConstant: public ValueDefinition {

public:

    /// Exemplar constructor.
    BoolConstant();

    virtual ~BoolConstant();

    /// Test if object can be replaced.
    //  Always false for constants.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual BoolConstant *clone(const std::string &name);

    /// Print the constant.
    virtual void print(std::ostream &) const;

    /// Print its value
    virtual void printValue(std::ostream &os) const;

    /// Return value.
    virtual bool getBool() const;

private:

    // Not implemented.
    BoolConstant(const BoolConstant &);
    void operator=(const BoolConstant &);

    // Clone constructor.
    BoolConstant(const std::string &name, BoolConstant *parent);
};

#endif // OPAL_BoolConstant_HH
