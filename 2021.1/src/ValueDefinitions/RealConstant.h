//
// Class RealConstant
//   The REAL CONSTANT definition.
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
#ifndef OPAL_RealConstant_HH
#define OPAL_RealConstant_HH

#include "AbstractObjects/ValueDefinition.h"


class RealConstant: public ValueDefinition {

public:

    /// Exemplar constructor.
    RealConstant();

    /// Constructor for built-in constants.
    RealConstant(const std::string &name, RealConstant *parent, double value);

    virtual ~RealConstant();

    /// Test if object can be replaced.
    //  Always false for constants.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual RealConstant *clone(const std::string &name);

    /// Print the constant.
    virtual void print(std::ostream &) const;

    /// Print its value
    virtual void printValue(std::ostream &os) const;

    /// Return value.
    virtual double getReal() const;

private:

    // Not implemented.
    RealConstant(const RealConstant &);
    void operator=(const RealConstant &);

    // Clone constructor.
    RealConstant(const std::string &name, RealConstant *parent);
};

#endif // OPAL_RealConstant_HH
