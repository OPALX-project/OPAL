//
// Class FM3DH5Block_nonscale
//   Class for dynamic non-scaled 3D field-maps stored in H5hut files.
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

#include "Fields/FM3DH5Block_nonscale.h"
#include "Physics/Physics.h"
#include "Physics/Units.h"
#include "Utilities/GeneralClassicException.h"

_FM3DH5Block_nonscale::_FM3DH5Block_nonscale (
    const std::string& filename
    ) : _Fieldmap (
        filename),
    _FM3DH5BlockBase (
        ) {
        Type = T3DDynamicH5Block;

        openFileMPIOCollective (filename);
        getFieldInfo ("Efield");
        getResonanceFrequency ();
        closeFile ();
}

_FM3DH5Block_nonscale::~_FM3DH5Block_nonscale (
    ) {
    freeMap();
}

FM3DH5Block_nonscale _FM3DH5Block_nonscale::create(const std::string& filename)
{
    return FM3DH5Block_nonscale(new _FM3DH5Block_nonscale(filename));
}

void _FM3DH5Block_nonscale::readMap (
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

    for (long long i = 0; i < num_gridpx_m * num_gridpy_m * num_gridpz_m; i++) {
        FieldstrengthEz_m[i] *= Units::MVpm2Vpm ;
        FieldstrengthEx_m[i] *= Units::MVpm2Vpm ;
        FieldstrengthEy_m[i] *= Units::MVpm2Vpm ;
        FieldstrengthHx_m[i] *= 1.0e6 * Physics::mu_0 ;
        FieldstrengthHy_m[i] *= 1.0e6 * Physics::mu_0 ;
        FieldstrengthHz_m[i] *= 1.0e6 * Physics::mu_0 ;
    }
    INFOMSG (level3
             << typeset_msg("3d dynamic (non-scaled) fieldmap '"
                            + Filename_m  + "' (H5hut format) read", "info")
             << endl);
}

void _FM3DH5Block_nonscale::freeMap (
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
}

bool _FM3DH5Block_nonscale::getFieldstrength (
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
        FieldstrengthHx_m, FieldstrengthHy_m, FieldstrengthHz_m, R);

    return false;
}
