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

#ifndef CLASSIC_FIELDMAP3DH5BLOCKBASE_H
#define CLASSIC_FIELDMAP3DH5BLOCKBASE_H

#include "Fields/Fieldmap.h"
#include <vector>

#include "H5hut.h"
static_assert (sizeof(double) == sizeof (h5_float64_t),
               "double and h5_float64_t are not the same type" );
static_assert (sizeof(long long) == sizeof (h5_int64_t),
               "long long and h5_int64_t are not the same type" );

class FM3DH5BlockBase: virtual public Fieldmap {

public:
    virtual void readMap (
        ) {};

    virtual void freeMap (
        ) {};

    virtual bool getFieldstrength (
        const Vector_t& /*R*/, Vector_t& /*E*/, Vector_t& /*B*/) const = 0;

    virtual void getFieldDimensions (
        double &zBegin, double &zEnd
        ) const {
        zBegin = zbegin_m;
        zEnd = zend_m;
    }

    virtual void getFieldDimensions (
        double &xIni, double &xFinal,
        double &yIni, double &yFinal,
        double &zIni, double &zFinal
        ) const {
        xIni = xbegin_m;
        xFinal = xend_m;
        yIni = ybegin_m;
        yFinal = yend_m;
        zIni = zbegin_m;
        zFinal = zend_m;
    }

    virtual bool getFieldDerivative (
        const Vector_t &/*R*/,
        Vector_t &/*E*/,
        Vector_t &/*B*/,
        const DiffDirection &/*dir*/
        ) const {
        return false;
    }

    virtual void swap(
        ) {};
    
    virtual void getInfo (
        Inform *msg);

    virtual double getFrequency (
        ) const;

    virtual void setFrequency (
        double freq);

    virtual void getOnaxisEz (
        std::vector<std::pair<double, double> >& F);

protected:
    FM3DH5BlockBase (
        ) {};

    virtual ~FM3DH5BlockBase (
        ) {};

    void openFileMPIOCollective (
        const std::string aFilename);

    long long getNumSteps (
        void);

    void setStep (
        const long long);

    void getFieldInfo (
        const char*);

    void getResonanceFrequency (
        void);

    void readField (
        const char* name,
        double* x,
        double* y,
        double* z
        );

    void closeFile (
        void);
    
    virtual bool isInside (
        const Vector_t &r
        ) const {
        return ((r(0) >= xbegin_m && r(0) < xend_m) &&
                (r(1) >= ybegin_m && r(1) < yend_m) &&
                (r(2) >= zbegin_m && r(2) < zend_m));
    }

    struct IndexTriplet {
        unsigned int i;
        unsigned int j;
        unsigned int k;
        Vector_t weight;
        IndexTriplet():
            i(0),
            j(0),
            k(0),
            weight(0.0)
        {}
    };

    /*
      The 3-dimensional fieldmaps are stored in a 1-dimensional arrays.
      Please note that the FORTRAN indexing scheme is used in H5hut!

      This functions maps the 3-dimensional index (i, j, k) to the
      corresponding index in the 1-dimensional array.
     */
    unsigned long getIndex (
        unsigned int i,
        unsigned int j,
        unsigned int k
        ) const {
        unsigned long result = j + k * num_gridpy_m;
        result = i + result * num_gridpx_m;
        
        return result;
    }

    IndexTriplet getIndex(const Vector_t &X) const {
        IndexTriplet idx;
        idx.i = std::floor((X(0) - xbegin_m) / hx_m);
        idx.j = std::floor((X(1) - ybegin_m) / hy_m);
        idx.k = std::floor((X(2) - zbegin_m) / hz_m);
        PAssert_LT(idx.i, num_gridpx_m - 1);
        PAssert_LT(idx.j, num_gridpy_m - 1);
        PAssert_LT(idx.k, num_gridpz_m - 1);

        idx.weight(0) = (X(0) - xbegin_m) / hx_m - idx.i;
        idx.weight(1) = (X(1) - ybegin_m) / hy_m - idx.j;
        idx.weight(2) = (X(2) - zbegin_m) / hz_m - idx.k;

        return idx;
    }

    double getWeightedData (
        const std::vector<double>& data,
        const IndexTriplet& idx,
        unsigned short corner) const;

    Vector_t interpolateTrilinearly (
        const std::vector<double>&,
        const std::vector<double>&,
        const std::vector<double>&,
        const Vector_t& X) const;

    enum : unsigned short {
        LX = 0,  // low X
        LY = 0,  // low Y
        LZ = 0,  // low Z
        HX = 4,  // high X
        HY = 2,  // high Y
        HZ = 1}; // high Z

    h5_file_t file_m;
    std::vector<double> FieldstrengthEz_m;    /**< 3D array with Ez */
    std::vector<double> FieldstrengthEx_m;    /**< 3D array with Ex */
    std::vector<double> FieldstrengthEy_m;    /**< 3D array with Ey */
    
    double xbegin_m;
    double xend_m;

    double ybegin_m;
    double yend_m;

    double zbegin_m;
    double zend_m;

    double hx_m;            /**< length between points in grid, x-direction */
    double hy_m;            /**< length between points in grid, y-direction */
    double hz_m;            /**< length between points in grid, z-direction */

    double num_gridpx_m;    /**< number of points after 0(not counted here) in grid, x-direction*/
    double num_gridpy_m;    /**< number of points after 0(not counted here) in grid, y-direction*/
    double num_gridpz_m;    /**< number of points after 0(not counted here) in grid, z-direction*/

    double frequency_m;

    bool swap_m;
    friend class Fieldmap;
};

#endif
