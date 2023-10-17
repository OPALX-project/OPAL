//
// Class OpalScalingFFAMagnet
//   The class provides the user interface for the SCALINGFFAMAGNET object.
//
// Copyright (c) 2017 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#ifndef OPAL_OPALSCALINGFFAMAGNET_H
#define OPAL_OPALSCALINGFFAMAGNET_H

#include "Elements/OpalBend.h"

/** OpalScalingFFAMagnet provides user interface information for the SCALINGFFA object
 *
 *  Defines three parameters - field map name, units for field, length for field
 */
class OpalScalingFFAMagnet: public OpalElement {

public:
    /** enum maps string to integer value for UI definitions */
    enum {
        B0 = COMMON,
        R0,
        FIELD_INDEX,
        TAN_DELTA,
        MAX_Y_POWER,
        END_FIELD_MODEL,
        END_LENGTH,
        CENTRE_LENGTH,
        RADIAL_NEG_EXTENT,
        RADIAL_POS_EXTENT,
        HEIGHT,
        MAGNET_START,
        MAGNET_END,
        AZIMUTHAL_EXTENT,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalScalingFFAMagnet();

    /** Destructor does nothing */
    virtual ~OpalScalingFFAMagnet();

    /** Inherited copy constructor */
    virtual OpalScalingFFAMagnet* clone(const std::string& name);

    /** Update the ScalingFFA with new parameters from UI parser */
    virtual void update();

private:
    // Not implemented.
    OpalScalingFFAMagnet(const OpalScalingFFAMagnet& );
    void operator=(const OpalScalingFFAMagnet& );

    // Clone constructor.
    OpalScalingFFAMagnet(const std::string& name, OpalScalingFFAMagnet* parent);

    void setupNamedEndField();
    void setupDefaultEndField();
};

#endif // OPAL_OPALSCALINGFFAMAGNET_H

