//
// Class FM3DH5Block
//   Class for dynamic 3D field-maps stored in H5hut files.
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

#ifndef CLASSIC_FIELDMAP3DH5BLOCK_H
#define CLASSIC_FIELDMAP3DH5BLOCK_H

#include "Fields/FM3DH5BlockBase.h"

#include <vector>

class _FM3DH5Block: public _FM3DH5BlockBase {

public:
    virtual bool getFieldstrength (
        const Vector_t &R, Vector_t &E, Vector_t &B) const;

    virtual ~_FM3DH5Block (
        );

private:
    _FM3DH5Block (
        const std::string& filename);

    static FM3DH5Block create(const std::string& filename);

    virtual void readMap (
        );

    virtual void freeMap (
        );

    std::vector<double> FieldstrengthHz_m;    /**< 3D array with Hz */
    std::vector<double> FieldstrengthHx_m;    /**< 3D array with Hx */
    std::vector<double> FieldstrengthHy_m;    /**< 3D array with Hy */

    friend class _Fieldmap;
    friend class _FM3DH5BlockBase;
};

using FM3DH5Block = std::shared_ptr<_FM3DH5Block>;

#endif