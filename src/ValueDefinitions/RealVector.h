//
// Class RealVector
//   The REAL VECTOR definition.
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
#ifndef OPAL_RealVector_HH
#define OPAL_RealVector_HH

#include "AbstractObjects/ValueDefinition.h"


class RealVector: public ValueDefinition {

public:

    /// Exemplar constructor.
    RealVector();

    virtual ~RealVector();

    /// Test for allowed replacement.
    //  True, if [b]rhs[/b] is a real vector.
    virtual bool canReplaceBy(Object *rhs);

    /// Make clone.
    virtual RealVector *clone(const std::string &name);

    /// Print the vector.
    virtual void print(std::ostream &) const;

    /// Print its value
    virtual void printValue(std::ostream &os) const;

    /// Return indexed value.
    virtual double getRealComponent(int) const;

private:

    // Not implemented.
    RealVector(const RealVector &);
    void operator=(const RealVector &);

    // Clone constructor.
    RealVector(const std::string &name, RealVector *parent);
};

#endif // OPAL_RealVector_HH
