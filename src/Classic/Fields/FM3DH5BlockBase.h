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
#include "Fields/FieldArray.h"

#include <cmath>
#include <vector>

#include "H5hut.h"
static_assert (sizeof(double) == sizeof (h5_float64_t),
               "double and h5_float64_t are not the same type" );
static_assert (sizeof(long long) == sizeof (h5_int64_t),
               "long long and h5_int64_t are not the same type" );

class _FM3DH5BlockBase: virtual public _Fieldmap {

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
    _FM3DH5BlockBase (
        ) {};

    virtual ~_FM3DH5BlockBase (
        ) {};

    void openFileMPIOCollective (
        const std::string& filename);

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

    /// Try to back three components with a read-only mapping of the file instead of
    /// reading them. All three succeed or none do. Returns false when the dataset layout
    /// is not directly mappable, in which case the caller reads as before.
    bool tryMapComponents (
        const char* name,
        long long step,
        FieldArray& x,
        FieldArray& y,
        FieldArray& z,
        std::size_t n
        );

    /// Collective over Ippl::getComm(): true if any sample of the three components is
    /// nonzero. Each rank scans an equal share, so every sample is read once per job
    /// and all ranks get the same answer.
    bool anyNonZero (
        const FieldArray& x,
        const FieldArray& y,
        const FieldArray& z
        ) const;

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

    /*
      Compute grid indices for a point X.

      Before calling this function a Isinside(X) test must
      to be performed!

      2-dim example with num_gridpx_m = 5, num_gridpy_m = 4

      3 +----+----+----+----+
        |    |    |    |    |
      2 +----+----+----+----+
        |    |    |    |   X|
      1 +----+----+----+----+
        |    |    |    |    |
      0 +----+----+----+----+
        0    1    2    3    4

      idx.x = 3
      idx.y = 1

      Notes:
      In above example: max for idx.x is num_gridpx_m - 2 which is 3.

      If X is close to a cell border, it can happen that the computation
      of the index is wrong due to rounding errors. In the example above
      idx.i == 4 instead of 3. To mitigate this issue we use __float128.
      To avoid out-of-bound errors we set the index to the max allowed
      value, if it is out of bound.

      std::fmod() doesn't support __float128 types! Best we can use is
      long double.
    */
    IndexTriplet getIndex(const Vector_t &X) const {
        IndexTriplet idx;

        // Derive the cell index and the in-cell weight from the SAME quantity, so the
        // two can never disagree, and keep the weight normalised to [0,1) as
        // getWeightedData() requires. Plain double throughout: the input X is already
        // double, so long double bought no accuracy, while fmodl() is a libm call and
        // x87 arithmetic does not vectorise.
        const double fx = (X(0) - xbegin_m) * inv_hx_m;
        const double fy = (X(1) - ybegin_m) * inv_hy_m;
        const double fz = (X(2) - zbegin_m) * inv_hz_m;

        // floor() before the cast: a negative fx would make the unsigned conversion
        // undefined. isInside() is supposed to exclude that, but this is one multiply
        // away from being free and removes the trap.
        const double ax = std::floor(fx);
        const double ay = std::floor(fy);
        const double az = std::floor(fz);

        idx.i = (unsigned int)std::min(std::max(ax, 0.0), (double)(num_gridpx_m - 2));
        idx.j = (unsigned int)std::min(std::max(ay, 0.0), (double)(num_gridpy_m - 2));
        idx.k = (unsigned int)std::min(std::max(az, 0.0), (double)(num_gridpz_m - 2));

        // Measured against idx.* rather than a*, so that when the clamp fires the
        // weight still describes the cell we actually selected.
        idx.weight(0) = fx - (double)idx.i;
        idx.weight(1) = fy - (double)idx.j;
        idx.weight(2) = fz - (double)idx.k;

        return idx;
    }

    double getWeightedData (
        const FieldArray& data,
        const IndexTriplet& idx,
        unsigned short corner) const;

    Vector_t interpolateTrilinearly (
        const FieldArray&,
        const FieldArray&,
        const FieldArray&,
        const Vector_t& X) const;

    enum : unsigned short {
        LX = 0,  // low X
        LY = 0,  // low Y
        LZ = 0,  // low Z
        HX = 4,  // high X
        HY = 2,  // high Y
        HZ = 1}; // high Z

    h5_file_t file_m;
    FieldArray FieldstrengthEz_m;    /**< 3D array with Ez */
    FieldArray FieldstrengthEx_m;    /**< 3D array with Ex */
    FieldArray FieldstrengthEy_m;    /**< 3D array with Ey */

    /// False once readMap() has found every Efield sample to be exactly zero; the E arrays
    /// are then released and never interpolated. Only readMap() implementations that scan
    /// clear it.
    bool hasE_m = true;

    double xbegin_m;
    double xend_m;

    double ybegin_m;
    double yend_m;

    double zbegin_m;
    double zend_m;

    double hx_m;            /**< length between points in grid, x-direction */
    double hy_m;            /**< length between points in grid, y-direction */
    double hz_m;            /**< length between points in grid, z-direction */

    // Reciprocals of the grid spacing, cached in getFieldInfo(). getIndex() runs
    // ~16x per particle per step; a division there is pure waste.
    double inv_hx_m;
    double inv_hy_m;
    double inv_hz_m;

    unsigned int num_gridpx_m;    /**< number of points after 0(not counted here) in grid, x-direction*/
    unsigned int num_gridpy_m;    /**< number of points after 0(not counted here) in grid, y-direction*/
    unsigned int num_gridpz_m;    /**< number of points after 0(not counted here) in grid, z-direction*/

    double frequency_m;

    bool swap_m;
    friend class _Fieldmap;
};

using FM3DH5BlockBase = std::shared_ptr<_FM3DH5BlockBase>;

#endif
