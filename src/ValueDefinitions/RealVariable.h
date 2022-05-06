//
// Class RealVariable
//   The REAL VARIABLE definition.
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
#ifndef OPAL_RealVariable_HH
#define OPAL_RealVariable_HH

#include "AbstractObjects/ValueDefinition.h"


class RealVariable: public ValueDefinition {

public:

    /// Exemplar constructor.
    RealVariable();

    /// Constructor for built-in variables.
    RealVariable(const std::string &name, RealVariable *parent, double value);

    virtual ~RealVariable();

    /// Test for allowed replacement.
    //  True, if [b]rhs[/b] is a real variable.
    virtual bool canReplaceBy(Object *rhs);

    /// Make clone.
    virtual RealVariable *clone(const std::string &name);

    /// Print the variable.
    virtual void print(std::ostream &) const;

    /// Print its value
    virtual void printValue(std::ostream &os) const;

    /// Return value.
    virtual double getReal() const;

private:

    // Not implemented.
    RealVariable(const RealVariable &);
    void operator=(const RealVariable &);

    // Clone constructor.
    RealVariable(const std::string &name, RealVariable *parent);
};

#endif // OPAL_RealVariable_HH
