//
// Class FM3DH5BlockBase
//   Base class for 3D field-maps in stored in H5hut files.
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

#include "Fields/FM3DH5BlockBase.h"
#include "Fields/Fieldmap.hpp"
#include "Physics/Physics.h"
#include "Utilities/GeneralClassicException.h"

void FM3DH5BlockBase::openFileMPIOCollective (
    const std::string aFilename
    ) {
    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    if (H5SetPropFileMPIOCollective (props, &comm) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::openFileMPIOCollective () ",
            "Cannot set MPIO collective!");
    }
    file_m = H5OpenFile (aFilename.c_str(), H5_O_RDONLY, props);
    if (file_m == (h5_file_t)H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::openFileMPIOCollective () ",
            "Cannot open file '" + aFilename + "'!");
    }
    H5CloseProp (props);
}

long long FM3DH5BlockBase::getNumSteps (void) {
    long long num_steps = H5GetNumSteps(file_m);
    if (num_steps <= 0) {
        if (num_steps == 0) {
            throw GeneralClassicException (
                "FM3DH5BlockBase::getNumSteps () ",
                "Number of time-steps in file '" + Filename_m + "' is zero!");
        } else {
            throw GeneralClassicException (
                "FM3DH5BlockBase::getNumSteps () ",
                "Query number of time-steps in file '" + Filename_m + "' failed!");
        }
    }
    return num_steps;
}

void FM3DH5BlockBase::setStep (long long step) {
    if (H5SetStep(file_m, step) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::setStep () ",
            "Cannot set time-step to " + std::to_string (step) +
            " in file '" + Filename_m + "'!");
    }
}
    
void FM3DH5BlockBase::getFieldInfo (const char* name) {
    long long last_step = getNumSteps () - 1;
    setStep (last_step);
    h5_size_t grid_rank;
    h5_size_t grid_dims[3];
    h5_size_t field_dims;
    h5_int64_t ftype;
    if (H5BlockGetFieldInfoByName (
                file_m, name,
                &grid_rank, grid_dims, &field_dims, &ftype
                ) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::GetFieldInfo () ",
            "Query of field info for " + std::string (name) +
            " in time-step " + std::to_string (last_step) +
            " in file '" + Filename_m + "' failed!");
    }
    num_gridpx_m = grid_dims[0];
    num_gridpy_m = grid_dims[1];
    num_gridpz_m = grid_dims[2];

    if (H5Block3dGetFieldSpacing (
            file_m, "Efield", &hx_m, &hy_m, &hz_m
            ) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::GetFieldInfo () ",
            "Query of field spacing"
            " in time-step " + std::to_string (last_step) +
            " in file '" + Filename_m + "' failed!");
    }
    
    if (H5Block3dGetFieldOrigin(
            file_m, "Efield", &xbegin_m, &ybegin_m, &zbegin_m) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::GetFieldInfo () ",
            "Query of field origin"
            " in time-step " + std::to_string (last_step) +
            " in file '" + Filename_m + "' failed!");
    }
    xend_m = xbegin_m + (num_gridpx_m - 1) * hx_m;
    yend_m = ybegin_m + (num_gridpy_m - 1) * hy_m;
    zend_m = zbegin_m + (num_gridpz_m - 1) * hz_m;
}

void FM3DH5BlockBase::getResonanceFrequency (void) {
    if (H5ReadFileAttribFloat64 (
            file_m, "Resonance Frequency(Hz)", &frequency_m
            ) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::GetResonanceFrequency () ",
            "Cannot read file attribute 'Resonance Frequency(Hz)'"
            " in file '" + Filename_m + "'!");
    }
    frequency_m *= Physics::two_pi;
}

void FM3DH5BlockBase::readField (
    const char* name,
    double* x,
    double* y,
    double* z
    ) {
    if (H5Block3dSetView(file_m,
                         0, num_gridpx_m - 1,
                         0, num_gridpy_m - 1,
                         0, num_gridpz_m - 1) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::ReadField () ",
            "Cannot set view "
            "0, " + std::to_string (num_gridpx_m) +
            "0, " + std::to_string (num_gridpy_m) +
            "0, " + std::to_string (num_gridpz_m) +
            " in file '" + Filename_m + "'!");
    }
    if (H5Block3dReadVector3dFieldFloat64 (
            file_m,
            name,
            x, y, z) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::ReadField () ",
            "Cannot read field " + std::string (name) +
            " in file '" + Filename_m + "'!");
    }
}

void FM3DH5BlockBase::closeFile (void) {
    if (H5CloseFile (file_m) == H5_ERR) {
        throw GeneralClassicException (
            "FM3DH5BlockBase::closeFile () ",
            "Error closing file '" + Filename_m + "'!");
    }
}

double FM3DH5BlockBase::getWeightedData (
    const std::vector<double>& data,
    const IndexTriplet& idx,
    unsigned short corner
    ) const {
    unsigned short switchX = ((corner & HX) >> 2);
    unsigned short switchY = ((corner & HY) >> 1);
    unsigned short switchZ =  (corner & HZ);
    double factorX = 0.5 + (1 - 2 * switchX) * (0.5 - idx.weight(0));
    double factorY = 0.5 + (1 - 2 * switchY) * (0.5 - idx.weight(1));
    double factorZ = 0.5 + (1 - 2 * switchZ) * (0.5 - idx.weight(2));

    unsigned long i = idx.i + switchX, j = idx.j + switchY, k = idx.k + switchZ;

    return factorX * factorY * factorZ * data[getIndex(i, j, k)];
}

Vector_t FM3DH5BlockBase::interpolateTrilinearly (
    const std::vector<double>& field_strength_x,
    const std::vector<double>& field_strength_y,
    const std::vector<double>& field_strength_z,
    const Vector_t& X
    ) const {
    IndexTriplet idx = getIndex (X);
    Vector_t result{0.0};

    result[0] = (getWeightedData(field_strength_x, idx, LX|LY|LZ) +
                 getWeightedData(field_strength_x, idx, LX|LY|HZ) +
                 getWeightedData(field_strength_x, idx, LX|HY|LZ) +
                 getWeightedData(field_strength_x, idx, LX|HY|HZ) +
                 getWeightedData(field_strength_x, idx, HX|LY|LZ) +
                 getWeightedData(field_strength_x, idx, HX|LY|HZ) +
                 getWeightedData(field_strength_x, idx, HX|HY|LZ) +
                 getWeightedData(field_strength_x, idx, HX|HY|HZ));

    result[1] = (getWeightedData(field_strength_y, idx, LX|LY|LZ) +
                 getWeightedData(field_strength_y, idx, LX|LY|HZ) +
                 getWeightedData(field_strength_y, idx, LX|HY|LZ) +
                 getWeightedData(field_strength_y, idx, LX|HY|HZ) +
                 getWeightedData(field_strength_y, idx, HX|LY|LZ) +
                 getWeightedData(field_strength_y, idx, HX|LY|HZ) +
                 getWeightedData(field_strength_y, idx, HX|HY|LZ) +
                 getWeightedData(field_strength_y, idx, HX|HY|HZ));

    result[2] = (getWeightedData(field_strength_z, idx, LX|LY|LZ) +
                 getWeightedData(field_strength_z, idx, LX|LY|HZ) +
                 getWeightedData(field_strength_z, idx, LX|HY|LZ) +
                 getWeightedData(field_strength_z, idx, LX|HY|HZ) +
                 getWeightedData(field_strength_z, idx, HX|LY|LZ) +
                 getWeightedData(field_strength_z, idx, HX|LY|HZ) +
                 getWeightedData(field_strength_z, idx, HX|HY|LZ) +
                 getWeightedData(field_strength_z, idx, HX|HY|HZ));

    return result;
}

void FM3DH5BlockBase::getInfo (Inform* msg) {
    (*msg) << Filename_m << " (3D dynamic) "
           << " xini= " << xbegin_m << " xfinal= " << xend_m
           << " yini= " << ybegin_m << " yfinal= " << yend_m
           << " zini= " << zbegin_m << " zfinal= " << zend_m << " [mm] " << endl;
    (*msg) << " hx= " << hx_m <<" hy= " << hy_m <<" hz= " << hz_m << " [mm] " <<endl;
}

double FM3DH5BlockBase::getFrequency () const {
    return frequency_m;
}

void FM3DH5BlockBase::setFrequency (double freq) {
    frequency_m = freq;
}

void FM3DH5BlockBase::getOnaxisEz (
    std::vector<std::pair<double, double>>& F
    ) {
    F.resize(num_gridpz_m);

    double Ez_max = 0.0;
    const double dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    const int index_x = -static_cast<int>(std::floor(xbegin_m / hx_m));
    const double lever_x = -xbegin_m / hx_m - index_x;

    const int index_y = -static_cast<int>(std::floor(ybegin_m / hy_m));
    const double lever_y = -ybegin_m / hy_m - index_y;
    long idx = index_x + index_y*num_gridpx_m + num_gridpx_m*num_gridpy_m;
    for (int i = 0;
         i < num_gridpz_m;
         i++, idx += num_gridpy_m * num_gridpx_m
        ) {
        F[i].first = dz * i;
        F[i].second =
            (1.0 - lever_x)   * (1.0 - lever_y) * FieldstrengthEz_m[idx]
            + lever_x         * (1.0 - lever_y) * FieldstrengthEz_m[idx + 1]
            + (1.0 - lever_x) * lever_y         * FieldstrengthEz_m[idx + num_gridpx_m]
            + lever_x         * lever_y         * FieldstrengthEz_m[idx + num_gridpx_m + 1];

        if(std::abs(F[i].second) > Ez_max) {
            Ez_max = std::abs(F[i].second);
        }
        F[i].second /= Ez_max;
    }
}
