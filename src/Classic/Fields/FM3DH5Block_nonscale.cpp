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

    // This class rescales the samples in place, so it always owns its buffers and is
    // never backed by a read-only mapping of the file.
    double* ex = FieldstrengthEx_m.mutableData();
    double* ey = FieldstrengthEy_m.mutableData();
    double* ez = FieldstrengthEz_m.mutableData();
    double* hx = FieldstrengthHx_m.mutableData();
    double* hy = FieldstrengthHy_m.mutableData();
    double* hz = FieldstrengthHz_m.mutableData();

    const double muFactor = 1.0e6 * Physics::mu_0;
    for (long long unsigned i = 0; i < num_gridpx_m * num_gridpy_m * num_gridpz_m; i++) {
        ez[i] *= Units::MVpm2Vpm ;
        ex[i] *= Units::MVpm2Vpm ;
        ey[i] *= Units::MVpm2Vpm ;
        hx[i] *= muFactor ;
        hy[i] *= muFactor ;
        hz[i] *= muFactor ;
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
    FieldstrengthEx_m.reset ();
    FieldstrengthEy_m.reset ();
    FieldstrengthEz_m.reset ();
    FieldstrengthHx_m.reset ();
    FieldstrengthHy_m.reset ();
    FieldstrengthHz_m.reset ();
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
