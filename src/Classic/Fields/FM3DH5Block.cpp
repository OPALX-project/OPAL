#include "Fields/FM3DH5Block.h"
#include "Fields/Fieldmap.hpp"
#include "H5hut.h"
#include "Physics/Physics.h"

#include <fstream>
#include <ios>

#include <assert.h>

using namespace std;
using Physics::mu_0;

FM3DH5Block::FM3DH5Block(std::string aFilename):Fieldmap(aFilename) {
    h5_err_t h5err;
#if defined (NDEBUG)
    (void)h5err;
#endif
    h5_size_t grid_rank;
    h5_size_t grid_dims[3];
    h5_size_t field_dims;
    char name[20];
    h5_size_t len_name = 20;
    h5_int64_t ftype;

    Type = T3DDynamicH5Block;

    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    h5err = H5SetPropFileMPIOCollective (props, &comm);
    assert (h5err != H5_ERR);
    h5_file_t file = H5OpenFile (aFilename.c_str(), H5_O_RDONLY, props);
    assert (file != (h5_file_t)H5_ERR);
    H5CloseProp (props);

    h5_int64_t last_step = H5GetNumSteps(file) - 1;
    assert (last_step >= 0);
    h5err = H5SetStep(file, last_step);
    assert (h5err != H5_ERR);

    h5_int64_t num_fields = H5BlockGetNumFields(file);
    assert (num_fields != H5_ERR);

    for(h5_ssize_t i = 0; i < num_fields; ++ i) {
        h5err = H5BlockGetFieldInfo(file, (h5_size_t)i, name, len_name, &grid_rank, grid_dims, &field_dims, &ftype);
        assert (h5err != H5_ERR);
        if(strcmp(name, "Efield") == 0) {
            num_gridpx_m = grid_dims[0];
            num_gridpy_m = grid_dims[1];
            num_gridpz_m = grid_dims[2];
        }
    }
    h5err = H5Block3dGetFieldSpacing(file, "Efield", &hx_m, &hy_m, &hz_m);
    assert (h5err != H5_ERR);
    h5err = H5Block3dGetFieldOrigin(file, "Efield", &xbegin_m, &ybegin_m, &zbegin_m);
    assert (h5err != H5_ERR);
    xend_m = xbegin_m + (num_gridpx_m - 1) * hx_m;
    yend_m = ybegin_m + (num_gridpy_m - 1) * hy_m;
    zend_m = zbegin_m + (num_gridpz_m - 1) * hz_m;

    //         xcentral_idx_m = static_cast<int>(fabs(xbegin_m) / hx_m);
    //         ycentral_idx_m = static_cast<int>(fabs(ybegin_m) / hy_m);


    h5err = H5ReadFileAttribFloat64(file, "Resonance Frequency(Hz)", &frequency_m);
    assert (h5err != H5_ERR);
    frequency_m *= Physics::two_pi;

    h5err = H5CloseFile(file);
    assert (h5err != H5_ERR);
}


FM3DH5Block::~FM3DH5Block() {
    freeMap();
}

void FM3DH5Block::readMap() {
    if (!FieldstrengthEz_m.empty()) {
        return;
    }
    h5_err_t h5err;
#if defined (NDEBUG)
    (void)h5err;
#endif
    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    h5err = H5SetPropFileMPIOCollective (props, &comm);
    assert (h5err != H5_ERR);
    h5_file_t file = H5OpenFile (Filename_m.c_str(), H5_O_RDONLY, props);
    assert (file != (h5_file_t)H5_ERR);
    H5CloseProp (props);

    long field_size = 0;

    h5_int64_t last_step = H5GetNumSteps(file) - 1;
    h5err = H5SetStep(file, last_step);
    assert (h5err != H5_ERR);

    field_size = (num_gridpx_m * num_gridpy_m * num_gridpz_m);
    FieldstrengthEx_m.resize(field_size);
    FieldstrengthEy_m.resize(field_size);
    FieldstrengthEz_m.resize(field_size);
    FieldstrengthHx_m.resize(field_size);
    FieldstrengthHy_m.resize(field_size);
    FieldstrengthHz_m.resize(field_size);
    h5err = H5Block3dSetView(file,
                             0, num_gridpx_m - 1,
                             0, num_gridpy_m - 1,
                             0, num_gridpz_m - 1);
    assert (h5err != H5_ERR);
    h5err = H5Block3dReadVector3dFieldFloat64(
        file,
        "Efield",
        &(FieldstrengthEx_m[0]),
        &(FieldstrengthEy_m[0]),
        &(FieldstrengthEz_m[0]));
    assert (h5err != H5_ERR);
    h5err = H5Block3dReadVector3dFieldFloat64(
        file,
        "Hfield",
        &(FieldstrengthHx_m[0]),
        &(FieldstrengthHy_m[0]),
        &(FieldstrengthHz_m[0]));
    assert (h5err != H5_ERR);

    h5err = H5CloseFile(file);
    assert (h5err != H5_ERR);

    INFOMSG(level3 << typeset_msg("fieldmap '" + Filename_m  + "' read", "info") << endl);
}

void FM3DH5Block::freeMap() {
    if(!FieldstrengthEz_m.empty()) {
        FieldstrengthEx_m.clear();
        FieldstrengthEy_m.clear();
        FieldstrengthEz_m.clear();
        FieldstrengthHx_m.clear();
        FieldstrengthHy_m.clear();
        FieldstrengthHz_m.clear();

        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << endl);
    }
}

bool FM3DH5Block::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    const int index_x = static_cast<int>(floor((R(0) - xbegin_m) / hx_m));
    const double lever_x = (R(0) - xbegin_m) / hx_m - index_x;

    const int index_y = static_cast<int>(floor((R(1) - ybegin_m) / hy_m));
    const double lever_y = (R(1) - ybegin_m) / hy_m - index_y;

    const int index_z = (int)floor((R(2) - zbegin_m) / hz_m);
    const double lever_z = (R(2) - zbegin_m) / hz_m - index_z;

    if((index_z < 0) || (index_z + 2 > num_gridpz_m)) {
        return false;
    }

    if(index_x < 0) {
        return false;
    }
    if(index_x + 2 > num_gridpx_m) {
        return false;
    }

    if(index_y < 0) {
        return false;
    }
    if(index_y + 2 > num_gridpy_m) {
        return false;
    }

    const long index1 = index_x + (index_y + index_z * num_gridpy_m) * num_gridpx_m;



    E(0) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEx_m[index1]
            + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEx_m[index1 + 1]
            + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEx_m[index1 + num_gridpx_m]
            + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEx_m[index1 + num_gridpx_m + 1]
            + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEx_m[index1 + num_gridpx_m * num_gridpy_m]
            + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEx_m[index1 + num_gridpx_m * num_gridpy_m + 1]
            + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
            + lever_x           * lever_y         * lever_z         * FieldstrengthEx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1];
    // msg<<R(0)<<", "<<R(1)<<", "<<R(2)<<", "<<index1<<", "<<index_x<<", "<<index_y<<", "<<index_z<<", "<<FieldstrengthEx_m[index1]<<", "<<FieldstrengthEy_m[index1]<<", "<<FieldstrengthEz_m[index1]<< endl;

    E(1) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEy_m[index1]
            + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEy_m[index1 + 1]
            + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEy_m[index1 + num_gridpx_m]
            + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEy_m[index1 + num_gridpx_m + 1]
            + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEy_m[index1 + num_gridpx_m * num_gridpy_m]
            + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEy_m[index1 + num_gridpx_m * num_gridpy_m + 1]
            + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
            + lever_x           * lever_y         * lever_z         * FieldstrengthEy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1];

    E(2) += (1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEz_m[index1]
            + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthEz_m[index1 + 1]
            + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthEz_m[index1 + num_gridpx_m]
            + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthEz_m[index1 + num_gridpx_m + 1]
            + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m]
            + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + 1]
            + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
            + lever_x           * lever_y         * lever_z         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1];

    B(0) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHx_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHx_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHx_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHx_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHx_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHx_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthHx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);

    B(1) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHy_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHy_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHy_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHy_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHy_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHy_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthHy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);

    B(2) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHz_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthHz_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthHz_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthHz_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthHz_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthHz_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthHz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthHz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);


    return false;
}

bool FM3DH5Block::getFieldDerivative(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {
    return false;
}

void FM3DH5Block::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = xbegin_m;
    rEnd = xend_m;
}
void FM3DH5Block::getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const {
    xIni = xbegin_m;
    xFinal = xend_m;
    yIni = ybegin_m;
    yFinal = yend_m;
    zIni = zbegin_m;
    zFinal = zend_m;
    //INFOMSG("xIni " << xIni << "xFinal " << xFinal << endl);
}
void FM3DH5Block::swap() { }

void FM3DH5Block::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (3D dynamic) "
           << " xini= " << xbegin_m << " xfinal= " << xend_m
           << " yini= " << ybegin_m << " yfinal= " << yend_m
           << " zini= " << zbegin_m << " zfinal= " << zend_m << " [mm] " << endl;
    (*msg) << " hx= " << hx_m <<" hy= " << hy_m <<" hz= " << hz_m << " [mm] " <<endl;
}

double FM3DH5Block::getFrequency() const {
    return frequency_m;
}

void FM3DH5Block::setFrequency(double freq) {
    frequency_m = freq;
}

void FM3DH5Block::getOnaxisEz(vector<pair<double, double> > & F) {
    double Ez_max = 0.0, dz = (zend_m - zbegin_m) / (num_gridpz_m - 1);
    const int index_x = -static_cast<int>(floor(xbegin_m / hx_m));
    const double lever_x = -xbegin_m / hx_m - index_x;

    const int index_y = -static_cast<int>(floor(ybegin_m / hy_m));
    const double lever_y = -ybegin_m / hy_m - index_y;

    long index1 = index_x + index_y * num_gridpx_m;

    F.resize(num_gridpz_m);

    for(int i = 0; i < num_gridpz_m; ++ i) {
        F[i].first = dz * i;
        F[i].second = (1.0 - lever_x)   * (1.0 - lever_y) * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m]
                      + lever_x           * (1.0 - lever_y) * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + 1]
                      + (1.0 - lever_x)   * lever_y         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
                      + lever_x           * lever_y         * FieldstrengthEz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1];


        if(fabs(F[i].second) > Ez_max) {
            Ez_max = fabs(F[i].second);
        }
        index1 += num_gridpy_m * num_gridpx_m;
    }

    for(int i = 0; i < num_gridpz_m; ++ i) {
        F[i].second /= Ez_max;
    }
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
