//
// Class OpalVerticalFFAMagnet
//   The class provides the user interface for the VERTICALFFA object
//
// Copyright (c) 2019 Chris Rogers
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
#ifndef OPAL_OPALVERTICALFFAMAGNET_H
#define OPAL_OPALVERTICALFFAMAGNET_H

#include "Elements/OpalElement.h"

class OpalVerticalFFAMagnet : public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        B0 = COMMON,
        FIELD_INDEX,
        WIDTH,
        MAX_HORIZONTAL_POWER,
        END_LENGTH,
        CENTRE_LENGTH,
        BB_LENGTH,
        HEIGHT_NEG_EXTENT,
        HEIGHT_POS_EXTENT,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalVerticalFFAMagnet();

    /** Destructor does nothing */
    virtual ~OpalVerticalFFAMagnet();

    /** Inherited copy constructor */
    virtual OpalVerticalFFAMagnet *clone(const std::string &name);

    /** Update the VerticalFFA with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalVerticalFFAMagnet(const OpalVerticalFFAMagnet &);
    void operator=(const OpalVerticalFFAMagnet &);

    // Clone constructor.
    OpalVerticalFFAMagnet(const std::string &name, OpalVerticalFFAMagnet *parent);

    static std::string docstring_m;
};

#endif // OPAL_OPALVERTICALFFAMAGNET_H

