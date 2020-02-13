#include "Fields/FM3DMagnetoStaticH5Block.h"
#include "Fields/Fieldmap.hpp"
#include "H5hut.h"
#include "Physics/Physics.h"

#include <fstream>
#include <ios>

#include <assert.h>

using namespace std;
using Physics::mu_0;

FM3DMagnetoStaticH5Block::FM3DMagnetoStaticH5Block(string aFilename):
    Fieldmap(aFilename) {
    h5_err_t h5err;
#if defined (NDEBUG)
    (void)h5err;
#endif
    h5_size_t grid_rank;
    h5_size_t grid_dims[3];
    h5_size_t field_dims;
    char name[20];
    h5_size_t len_name = 20;

    Type = T3DMagnetoStaticH5Block;

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
        h5_int64_t ftype;
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


FM3DMagnetoStaticH5Block::~FM3DMagnetoStaticH5Block() {
    freeMap();
}

void FM3DMagnetoStaticH5Block::readMap() {
    if (!FieldstrengthEz_m.empty()) {
        return;
    }
    h5_int64_t h5err;
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
    int Nnodes = Ippl::getNodes();//min(20, Ippl::getNodes());
    int Nz_avrg = static_cast<int>(floor(0.5 + num_gridpz_m / Nnodes));
    int Nz_diff = Nz_avrg * Nnodes - num_gridpz_m;
    int signNz = Nz_diff > 0 ? 1 : -1;
    int *Nz_read_start = new int[Ippl::getNodes() + 1];
    int *Nz_read_length = new int[Ippl::getNodes()];
    int N_read_start;
    int start = 0;
    // int rbuf_size;

    h5_int64_t last_step = H5GetNumSteps(file) - 1;
    assert (last_step >= 0);
    h5err = H5SetStep(file, last_step);
    assert (h5err != H5_ERR);

    for(int i = 0; i < abs(Nz_diff); ++ i) {
        Nz_read_length[i] = Nz_avrg - signNz;
        Nz_read_start[i] = start;
        start += Nz_read_length[i];
    }
    for(int i = abs(Nz_diff); i < Nnodes; ++ i) {
        Nz_read_length[i] = Nz_avrg;
        Nz_read_start[i] = start;
        start += Nz_read_length[i];
    }
    for(int i = Nnodes; i < Ippl::getNodes(); ++ i) {
        Nz_read_length[i] = 0;
        Nz_read_start[i] = start;
    }
    Nz_read_start[Ippl::getNodes()] = start;

    N_read_start = Nz_read_start[Ippl::myNode()] * num_gridpx_m * num_gridpy_m;

    // rbuf_size = max(Nz_avrg, Nz_avrg - signNz);
    // std::unique_ptr<double> rbuf(new double[Ippl::getNodes() * rbuf_size]);

    h5err = H5Block3dSetView(file,
                             0, num_gridpx_m - 1,
                             0, num_gridpy_m - 1,
                             Nz_read_start[Ippl::myNode()], Nz_read_start[Ippl::myNode() + 1] - 1);
    assert (h5err != H5_ERR);

    field_size = (num_gridpx_m * num_gridpy_m * num_gridpz_m);
    FieldstrengthEx_m.resize(field_size);
    FieldstrengthEy_m.resize(field_size);
    FieldstrengthEz_m.resize(field_size);
    FieldstrengthBx_m.resize(field_size);
    FieldstrengthBy_m.resize(field_size);
    FieldstrengthBz_m.resize(field_size);

    h5err = H5Block3dReadVector3dFieldFloat64(
        file,
        "Efield",
        &(FieldstrengthEx_m[N_read_start]),
        &(FieldstrengthEy_m[N_read_start]),
        &(FieldstrengthEz_m[N_read_start]));
    assert (h5err != H5_ERR);

    h5err = H5Block3dReadVector3dFieldFloat64(
        file,
        "Bfield",
        &(FieldstrengthBx_m[N_read_start]),
        &(FieldstrengthBy_m[N_read_start]),
        &(FieldstrengthBz_m[N_read_start]));
    assert (h5err != H5_ERR);

    for(int i = 0; i < Nnodes; ++ i) {
        int N_read_start = Nz_read_start[i] * num_gridpx_m * num_gridpy_m;
        int N_read_length = Nz_read_length[i] * num_gridpx_m * num_gridpy_m;
        MPI_Bcast(&(FieldstrengthEx_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
        MPI_Bcast(&(FieldstrengthEy_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
        MPI_Bcast(&(FieldstrengthEz_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
        MPI_Bcast(&(FieldstrengthBx_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
        MPI_Bcast(&(FieldstrengthBy_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
        MPI_Bcast(&(FieldstrengthBz_m[N_read_start]), N_read_length, MPI_DOUBLE, i, Ippl::getComm());
    }

    h5err = H5CloseFile(file);
    assert (h5err != H5_ERR);

    delete[] Nz_read_start;
    delete[] Nz_read_length;

    INFOMSG(level3 << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl);

}

void FM3DMagnetoStaticH5Block::freeMap() {
    if(!FieldstrengthEz_m.empty()) {
        FieldstrengthEx_m.clear();
        FieldstrengthEy_m.clear();
        FieldstrengthEz_m.clear();
        FieldstrengthBx_m.clear();
        FieldstrengthBy_m.clear();
        FieldstrengthBz_m.clear();
        INFOMSG(level3 << typeset_msg("freed fieldmap '" + Filename_m + "'", "info") << "\n"
                << endl);
    }

}

bool FM3DMagnetoStaticH5Block::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    const int index_x = static_cast<int>(floor((R(0) - xbegin_m) / hx_m));
    const double lever_x = (R(0) - xbegin_m) / hx_m - index_x;

    const int index_y = static_cast<int>(floor((R(1) - ybegin_m) / hy_m));
    const double lever_y = (R(1) - ybegin_m) / hy_m - index_y;

    const int index_z = (int)floor((R(2) - zbegin_m) / hz_m);
    const double lever_z = (R(2) - zbegin_m) / hz_m - index_z;

    if((index_z < 0) || (index_z + 2 > num_gridpz_m)) {
        return false;
    }

    if(index_x + 2 > num_gridpx_m || index_y + 2 > num_gridpy_m) {
        ERRORMSG("Field strength at R " << R << " out of field maps" << endl);
        return false;// if x or y > xmax or ymax, do nothing just let E and B unchange in the element with this field map, we check if a particle outside domain in BoundaryGeometry, not here. That's resonable if we have field mapped without covering entire boundary of an element.
    }
    if(index_x < 0 || index_y < 0) {
        ERRORMSG("Field strength at R " << R << " out of field maps" << endl);
        return false;// if x or y < xmin or ymin, do nothing just let E and B unchange in the element with this field map, we check if a particle outside domain in BoundaryGeometry, not here. That's resonable if we have field mapped without covering entire boundary of an element.
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

    B(0) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthBx_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthBx_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthBx_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthBx_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthBx_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthBx_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthBx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthBx_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);

    B(1) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthBy_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthBy_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthBy_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthBy_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthBy_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthBy_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthBy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthBy_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);

    B(2) += ((1.0 - lever_x) * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthBz_m[index1]
             + lever_x           * (1.0 - lever_y) * (1.0 - lever_z) * FieldstrengthBz_m[index1 + 1]
             + (1.0 - lever_x)   * lever_y         * (1.0 - lever_z) * FieldstrengthBz_m[index1 + num_gridpx_m]
             + lever_x           * lever_y         * (1.0 - lever_z) * FieldstrengthBz_m[index1 + num_gridpx_m + 1]
             + (1.0 - lever_x)   * (1.0 - lever_y) * lever_z         * FieldstrengthBz_m[index1 + num_gridpx_m * num_gridpy_m]
             + lever_x           * (1.0 - lever_y) * lever_z         * FieldstrengthBz_m[index1 + num_gridpx_m * num_gridpy_m + 1]
             + (1.0 - lever_x)   * lever_y         * lever_z         * FieldstrengthBz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m]
             + lever_x           * lever_y         * lever_z         * FieldstrengthBz_m[index1 + num_gridpx_m * num_gridpy_m + num_gridpx_m + 1]);


    return false;

}

bool FM3DMagnetoStaticH5Block::getFieldDerivative(const Vector_t &/*R*/, Vector_t &/*E*/, Vector_t &/*B*/, const DiffDirection &/*dir*/) const {
    return false;
}

void FM3DMagnetoStaticH5Block::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = xbegin_m;
    rEnd = xend_m;
}
void FM3DMagnetoStaticH5Block::getFieldDimensions(double &/*xIni*/, double &/*xFinal*/, double &/*yIni*/, double &/*yFinal*/, double &/*zIni*/, double &/*zFinal*/) const {}

void FM3DMagnetoStaticH5Block::swap() {}

void FM3DMagnetoStaticH5Block::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (3DMagnetoStaticH5Block); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM3DMagnetoStaticH5Block::getFrequency() const {
    return 0.0;
}

void FM3DMagnetoStaticH5Block::setFrequency(double /*freq*/)
{ ;}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode:nil
// End:
