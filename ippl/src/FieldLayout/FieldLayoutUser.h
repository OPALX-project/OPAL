//
// FieldLayoutUser is a base class for all classes which need to use
// a FieldLayout - it is derived from User, which provides a virtual
// function 'notifyUserOfDelete' which is called when the FieldLayout
// is deleted, and the virtual function 'Repartition' which is called
// when a Field needs to be redistributed between processors.
//
// Copyright (c) 2003 - 2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
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

#ifndef FIELD_LAYOUT_USER_H
#define FIELD_LAYOUT_USER_H

#include "Utility/User.h"
#include "Utility/UserList.h"

// class definition
class FieldLayoutUser : public User {
public:
    // constructor - the base class selects a unique ID value
    FieldLayoutUser(){};

    // destructor, nothing to do here
    virtual ~FieldLayoutUser(){};

    //
    // virtual functions for FieldLayoutUser's
    //

    // Repartition onto a new layout
    virtual void Repartition(UserList *) = 0;
};

#endif  // FIELD_LAYOUT_USER_H
