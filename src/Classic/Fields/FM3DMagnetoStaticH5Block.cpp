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

FM3DMagnetoStaticH5Block::FM3DMagnetoStaticH5Block (
    std::string aFilename
    ) : Fieldmap (
        aFilename),
    FM3DH5BlockBase (
        ) {
        Type = T3DMagnetoStaticH5Block;

        openFileMPIOCollective (aFilename);
        getFieldInfo ("Efield");
        getResonanceFrequency ();
        closeFile ();
}

FM3DMagnetoStaticH5Block::~FM3DMagnetoStaticH5Block (
    ) {
    freeMap ();
}

void FM3DMagnetoStaticH5Block::readMap (
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
    FieldstrengthBx_m.resize (field_size);
    FieldstrengthBy_m.resize (field_size);
    FieldstrengthBz_m.resize (field_size);

    readField (
        "Efield",
        &(FieldstrengthEx_m[0]),
        &(FieldstrengthEy_m[0]),
        &(FieldstrengthEz_m[0]));
    readField (
        "Bfield",
        &(FieldstrengthBx_m[0]),
        &(FieldstrengthBy_m[0]),
        &(FieldstrengthBz_m[0]));

    closeFile ();
    INFOMSG (level3
             << typeset_msg("3d magneto static fieldmap '"
                            + Filename_m  + "' (H5hut format) read", "info")
             << endl);
}
void FM3DMagnetoStaticH5Block::freeMap (
    ) {
    if(FieldstrengthEz_m.empty()) {
        return;
    }
    FieldstrengthEx_m.clear();
    FieldstrengthEy_m.clear();
    FieldstrengthEz_m.clear();
    FieldstrengthBx_m.clear();
    FieldstrengthBy_m.clear();
    FieldstrengthBz_m.clear();
    INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m + "'", "info")
            << endl);
}

bool FM3DMagnetoStaticH5Block::getFieldstrength (
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

double FM3DMagnetoStaticH5Block::getFrequency (
    ) const {
    return 0.0;
}
