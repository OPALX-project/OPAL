//
// Class _FM3DH5Block
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

#include "Fields/FM3DH5Block.h"

_FM3DH5Block::_FM3DH5Block (
    const std::string& filename
    ) : _Fieldmap(
        filename
        ) {
        Type = T3DDynamicH5Block;

        openFileMPIOCollective (filename);
        getFieldInfo ("Efield");
        getResonanceFrequency ();
        closeFile ();
}

_FM3DH5Block::~_FM3DH5Block (
    ) {
    freeMap ();
}

FM3DH5Block _FM3DH5Block::create(const std::string& filename)
{
    return FM3DH5Block(new _FM3DH5Block(filename));
}

void _FM3DH5Block::readMap (
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
        if (tryMapComponents ("Hfield", last_step,
                              FieldstrengthHx_m, FieldstrengthHy_m,
                              FieldstrengthHz_m, field_size)) {
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
    FieldstrengthHx_m.allocate (field_size);
    FieldstrengthHy_m.allocate (field_size);
    FieldstrengthHz_m.allocate (field_size);

    readField (
        "Efield",
        FieldstrengthEx_m.mutableData(),
        FieldstrengthEy_m.mutableData(),
        FieldstrengthEz_m.mutableData());
    readField (
        "Hfield",
        FieldstrengthHx_m.mutableData(),
        FieldstrengthHy_m.mutableData(),
        FieldstrengthHz_m.mutableData());

    closeFile ();
    INFOMSG (level3
             << typeset_msg("3d dynamic fieldmap '"
                            + Filename_m  + "' (H5hut format) read", "info")
             << endl);
}

void _FM3DH5Block::freeMap (
    ) {
    if(FieldstrengthEz_m.empty ()) {
        return;
    }
    FieldstrengthEx_m.reset ();
    FieldstrengthEy_m.reset ();
    FieldstrengthEz_m.reset ();
    FieldstrengthHx_m.reset ();
    FieldstrengthHy_m.reset ();
    FieldstrengthHz_m.reset ();
}

bool _FM3DH5Block::getFieldstrength (
    const Vector_t& R,
    Vector_t& E,
    Vector_t& B
    ) const {
    if (!isInside(R)) {
        return true;
    }
    E += interpolateTrilinearly (FieldstrengthEx_m, FieldstrengthEy_m, FieldstrengthEz_m, R);
    B += interpolateTrilinearly (FieldstrengthHx_m, FieldstrengthHy_m, FieldstrengthHz_m, R);
    return false;
}
