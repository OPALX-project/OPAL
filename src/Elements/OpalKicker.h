//
// Class OpalKicker
//   The KICKER element.
//   Note the sign convention:  Positive kicks bend particles to positive x or
//   y respectively.
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
#ifndef OPAL_OpalKicker_HH
#define OPAL_OpalKicker_HH

#include "Elements/OpalElement.h"


class OpalKicker: public OpalElement {

public:

    /// The attributes of class OpalKicker.
    enum {
        HKICK = COMMON,  // The horizontal kicker strength.
        VKICK,           // The vertical kicker strength.
        DESIGNENERGY,    // The mean kinetic energy at exit
        K0,              // The normal dipole field
        K0S,             // The skew dipole field
        SIZE
    };

    /// Exemplar constructor.
    OpalKicker();

    virtual ~OpalKicker();

    /// Make clone.
    virtual OpalKicker *clone(const std::string &name);


    // JMJ 18/12/2000 Following method not needed, commented out, delete after next CVS commit.
    //BEGIN JMJ 15/12/2000, adding missing print method
    // Print the kicker
    //  Handle printing in OPAL-8 format.
    //  virtual void print(std::ostream &) const;
    //END   JMJ 15/12/2000, adding missing print method

    /// Update the embedded CLASSIC corrector.
    virtual void update();

private:

    // Not implemented.
    OpalKicker(const OpalKicker &);
    void operator=(const OpalKicker &);

    // Clone constructor.
    OpalKicker(const std::string &name, OpalKicker *parent);
};

#endif // OPAL_OpalKicker_HH
