//
// Class AttList
//   The ATTLIST command.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef OPAL_AttList_HH
#define OPAL_AttList_HH

#include "AbstractObjects/Action.h"
#include <iosfwd>

class Beamline;


class AttList: public Action {

public:

    /// Exemplar constructor.
    AttList();

    virtual ~AttList();

    /// Make clone.
    virtual AttList *clone(const std::string &name);

    /// Execute the command.
    virtual void execute();

private:

    // Not implemented.
    AttList(const AttList &);
    void operator=(const AttList &);

    // Clone constructor.
    AttList(const std::string &name, AttList *parent);

    // The working routine.
    void writeTable(const Beamline &line, std::ostream &os);
};

#endif // OPAL_AttList_HH
