//
// Class FM3DMagnetoStaticH5Block
//   Class for magneto-static 3D field-maps stored in H5hut files.
//
// Copyright (c) 2020, Achim Gsell, Paul Scherrer Institut, Villigen PSI, Switzerland
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
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//

#ifndef CLASSIC_FIELDMAP3DMAGNETOSTATICH5BLOCK_H
#define CLASSIC_FIELDMAP3DMAGNETOSTATICH5BLOCK_H

#include "Fields/FM3DH5BlockBase.h"

#include <vector>

class FM3DMagnetoStaticH5Block: public FM3DH5BlockBase {

public:
    virtual bool getFieldstrength (
        const Vector_t &R, Vector_t &E, Vector_t &B) const;
    
private:
    FM3DMagnetoStaticH5Block (
        std::string aFilename);

    virtual ~FM3DMagnetoStaticH5Block (
        );

    virtual void readMap (
        );

    virtual void freeMap (
        );

    virtual double getFrequency (
        ) const;

    std::vector<double> FieldstrengthBz_m;    /**< 3D array with Bz */
    std::vector<double> FieldstrengthBx_m;    /**< 3D array with Bx */
    std::vector<double> FieldstrengthBy_m;    /**< 3D array with By */

    friend class Fieldmap;
    friend class FM3DH5BlockBase;
};

#endif
