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

#include "Fields/FM3DMagnetoStaticH5Block.h"
#include "Utilities/GeneralClassicException.h"

_FM3DMagnetoStaticH5Block::_FM3DMagnetoStaticH5Block (const std::string& filename)
    : _Fieldmap(filename)
    , _FM3DH5BlockBase () {
    Type = T3DMagnetoStaticH5Block;

    openFileMPIOCollective (filename);
    getFieldInfo ("Efield");
    getResonanceFrequency ();
    closeFile ();
}

_FM3DMagnetoStaticH5Block::~_FM3DMagnetoStaticH5Block () {
    freeMap ();
}

FM3DMagnetoStaticH5Block _FM3DMagnetoStaticH5Block::create(const std::string& filename)
{
    return FM3DMagnetoStaticH5Block(new _FM3DMagnetoStaticH5Block(filename));
}

void _FM3DMagnetoStaticH5Block::readMap (
    ) {
    if (!FieldstrengthEz_m.empty()) {
        return;
    }
    openFileMPIOCollective (Filename_m);
    long long last_step = getNumSteps () - 1;
    setStep (last_step);

    size_t field_size = num_gridpx_m * num_gridpy_m * num_gridpz_m;

    // One physical copy per node instead of one per rank, when the file allows it.
    if (tryMapComponents ("Efield", last_step,
                          FieldstrengthEx_m, FieldstrengthEy_m, FieldstrengthEz_m, field_size)) {
        if (tryMapComponents ("Bfield", last_step,
                              FieldstrengthBx_m, FieldstrengthBy_m,
                              FieldstrengthBz_m, field_size)) {
            closeFile ();
            return;
        }
        FieldstrengthEx_m.reset ();
        FieldstrengthEy_m.reset ();
        FieldstrengthEz_m.reset ();
    }
    FieldstrengthEx_m.allocate (field_size);
    FieldstrengthEy_m.allocate (field_size);
    FieldstrengthEz_m.allocate (field_size);
    FieldstrengthBx_m.allocate (field_size);
    FieldstrengthBy_m.allocate (field_size);
    FieldstrengthBz_m.allocate (field_size);

    readField (
        "Efield",
        FieldstrengthEx_m.mutableData(),
        FieldstrengthEy_m.mutableData(),
        FieldstrengthEz_m.mutableData());
    readField (
        "Bfield",
        FieldstrengthBx_m.mutableData(),
        FieldstrengthBy_m.mutableData(),
        FieldstrengthBz_m.mutableData());

    closeFile ();
    INFOMSG (level3
             << typeset_msg("3d magneto static fieldmap '"
                            + Filename_m  + "' (H5hut format) read", "info")
             << endl);
}
void _FM3DMagnetoStaticH5Block::freeMap (
    ) {
    if(FieldstrengthEz_m.empty()) {
        return;
    }
    FieldstrengthEx_m.reset ();
    FieldstrengthEy_m.reset ();
    FieldstrengthEz_m.reset ();
    FieldstrengthBx_m.reset ();
    FieldstrengthBy_m.reset ();
    FieldstrengthBz_m.reset ();
}

bool _FM3DMagnetoStaticH5Block::getFieldstrength (
    const Vector_t& R,
    Vector_t& E,
    Vector_t& B
    ) const {
    if (!isInside(R)) {
        return true;
    }
    E += interpolateTrilinearly (
        FieldstrengthEx_m, FieldstrengthEy_m, FieldstrengthEz_m, R);
    B += interpolateTrilinearly (
        FieldstrengthBx_m, FieldstrengthBy_m, FieldstrengthBz_m, R);

    return false;
}

double _FM3DMagnetoStaticH5Block::getFrequency (
    ) const {
    return 0.0;
}
