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

#include "Fields/FM3DH5Block.h"

FM3DH5Block::FM3DH5Block (
    std::string aFilename
    ) : Fieldmap (
        aFilename
        ) {
        Type = T3DDynamicH5Block;

        openFileMPIOCollective (aFilename);
        getFieldInfo ("Efield");
        getResonanceFrequency ();
        closeFile ();
}

FM3DH5Block::~FM3DH5Block (
    ) {
    freeMap ();
}

void FM3DH5Block::readMap (
    ) {
    if (!FieldstrengthEz_m.empty()) {
        return;
    }
    openFileMPIOCollective (Filename_m);
    long long last_step = getNumSteps () - 1;
    setStep (last_step);

    size_t field_size = num_gridpx_m * num_gridpy_m * num_gridpz_m;
    FieldstrengthEx_m.resize (field_size);
    FieldstrengthEy_m.resize (field_size);
    FieldstrengthEz_m.resize (field_size);
    FieldstrengthHx_m.resize (field_size);
    FieldstrengthHy_m.resize (field_size);
    FieldstrengthHz_m.resize (field_size);

    readField (
        "Efield",
        &(FieldstrengthEx_m[0]),
        &(FieldstrengthEy_m[0]),
        &(FieldstrengthEz_m[0]));
    readField (
        "Hfield",
        &(FieldstrengthHx_m[0]),
        &(FieldstrengthHy_m[0]),
        &(FieldstrengthHz_m[0]));

    closeFile ();
    INFOMSG (level3
             << typeset_msg("3d dynamic fieldmap '"
                            + Filename_m  + "' (H5hut format) read", "info")
             << endl);
}

void FM3DH5Block::freeMap (
    ) {
    if(FieldstrengthEz_m.empty ()) {
        return;
    }
    FieldstrengthEx_m.clear ();
    FieldstrengthEy_m.clear ();
    FieldstrengthEz_m.clear ();
    FieldstrengthHx_m.clear ();
    FieldstrengthHy_m.clear ();
    FieldstrengthHz_m.clear ();

    INFOMSG (level3
             << typeset_msg ("freed fieldmap '" + Filename_m + "'", "info")
             << endl);
}

bool FM3DH5Block::getFieldstrength (
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
