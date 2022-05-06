//
// Declaration of the BoundaryGeometry class
//
// Copyright (c) 200x - 2020, Achim Gsell,
//                            Paul Scherrer Institut, Villigen PSI, Switzerland
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
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

//#define   ENABLE_DEBUG

#include "Structure/BoundaryGeometry.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <string>

#include "H5hut.h"
#include <cfloat>

#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunchBase.h"
#include "Expressions/SRefExpr.h"
#include "Elements/OpalBeamline.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"

#include <boost/filesystem.hpp>

#include <gsl/gsl_sys.h>

extern Inform* gmsg;

#define SQR(x) ((x)*(x))
#define PointID(triangle_id, vertex_id) Triangles_m[triangle_id][vertex_id]
#define Point(triangle_id, vertex_id)   Points_m[Triangles_m[triangle_id][vertex_id]]

/*
  In the following namespaces various approximately floating point
  comparisons are implemented. The used implementation is selected
  via

  namespaces cmp = IMPLEMENTATION;

*/

/*
  First we define some macros for function common in all namespaces.
*/
#define FUNC_EQ(x, y) inline bool eq(double x, double y) { \
        return almost_eq(x, y);                            \
    }

#define FUNC_EQ_ZERO(x) inline bool eq_zero(double x) { \
        return almost_eq_zero(x);                       \
    }

#define FUNC_LE(x, y) inline bool le(double x, double y) { \
        if (almost_eq(x, y)) {                             \
            return true;                                   \
        }                                                  \
        return x < y;                                      \
    }

#define FUNC_LE_ZERO(x) inline bool le_zero(double x) { \
        if (almost_eq_zero(x)) {                        \
            return true;                                \
        }                                               \
        return x < 0.0;                                 \
    }

#define FUNC_LT(x, y) inline bool lt(double x, double y) { \
        if (almost_eq(x, y)) {                             \
            return false;                                  \
        }                                                  \
        return x < y;                                      \
    }

#define FUNC_LT_ZERO(x) inline bool lt_zero(double x) {      \
        if (almost_eq_zero(x)) {                             \
            return false;                                    \
        }                                                    \
        return x < 0.0;                                      \
    }

#define FUNC_GE(x, y) inline bool ge(double x, double y) {      \
        if (almost_eq(x, y)) {                                  \
            return true;                                        \
        }                                                       \
        return x > y;                                           \
    }

#define FUNC_GE_ZERO(x) inline bool ge_zero(double x) { \
        if (almost_eq_zero(x)) {                        \
            return true;                                \
        }                                               \
        return x > 0.0;                                 \
    }

#define FUNC_GT(x, y) inline bool gt(double x, double y) { \
        if (almost_eq(x, y)) {                             \
            return false;                                  \
        }                                                  \
        return x > y;                                      \
    }

#define FUNC_GT_ZERO(x) inline bool gt_zero(double x) { \
        if (almost_eq_zero(x)) {                        \
            return false;                               \
        }                                               \
        return x > 0.0;                                 \
    }

namespace cmp_diff {

    /*
      Link:
      https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    */
    inline bool almost_eq(double A, double B, double maxDiff = 1e-15, double maxRelDiff = DBL_EPSILON) {
        // Check if the numbers are really close -- needed
        // when comparing numbers near zero.
        const double diff = std::abs(A - B);
        if (diff <= maxDiff)
            return true;

        A = std::abs(A);
        B = std::abs(B);
        const double largest = (B > A) ? B : A;

        if (diff <= largest * maxRelDiff)
            return true;
        return false;
    }

    inline bool almost_eq_zero(double A, double maxDiff = 1e-15) {
        const double diff = std::abs(A);
        return (diff <= maxDiff);
    }

    FUNC_EQ(x, y);
    FUNC_EQ_ZERO(x);
    FUNC_LE(x, y);
    FUNC_LE_ZERO(x);
    FUNC_LT(x, y);
    FUNC_LT_ZERO(x);
    FUNC_GE(x, y);
    FUNC_GE_ZERO(x);
    FUNC_GT(x, y);
    FUNC_GT_ZERO(x);
}

namespace cmp_ulp_obsolete {
    /*
      See:
      https://www.cygnus-software.com/papers/comparingfloats/comparing_floating_point_numbers_obsolete.htm
    */
    inline bool almost_eq(double A, double B, double maxDiff = 1e-20, int maxUlps = 1000) {
        // Make sure maxUlps is non-negative and small enough that the
        // default NAN won't compare as equal to anything.
        // assert(maxUlps > 0 && maxUlps < 4 * 1024 * 1024);

        // handle NaN's
        // Note: comparing something with a NaN is always false!
        if (std::isnan(A) || std::isnan(B)) {
            return false;
        }

        if (std::abs (A - B) <= maxDiff) {
            return true;
        }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
        auto aInt = *(int64_t*)&A;
#pragma GCC diagnostic pop
        // Make aInt lexicographically ordered as a twos-complement int
        if (aInt < 0) {
            aInt = 0x8000000000000000 - aInt;
        }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
        auto bInt = *(int64_t*)&B;
#pragma GCC diagnostic pop
        // Make bInt lexicographically ordered as a twos-complement int
        if (bInt < 0) {
            bInt = 0x8000000000000000 - bInt;
        }

        if (std::abs (aInt - bInt) <= maxUlps) {
            return true;
        }
        return false;
    }

    inline bool almost_eq_zero(double A, double maxDiff = 1e-15) {
        // no need to handle NaN's!
        return (std::abs(A) <= maxDiff);
    }
    FUNC_EQ(x, y);
    FUNC_EQ_ZERO(x);
    FUNC_LE(x, y);
    FUNC_LE_ZERO(x);
    FUNC_LT(x, y);
    FUNC_LT_ZERO(x);
    FUNC_GE(x, y);
    FUNC_GE_ZERO(x);
    FUNC_GT(x, y);
    FUNC_GT_ZERO(x);
}

namespace cmp_ulp {
    /*
      See:
      https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
    */


    inline bool almost_eq (double A, double B, double maxDiff = 1e-20, int maxUlps = 1000) {
        // handle NaN's
        if (std::isnan (A) || std::isnan (B)) {
            return false;
        }

        // Check if the numbers are really close -- needed
        // when comparing numbers near zero.
        if (std::abs (A - B) <= maxDiff)
            return true;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
        auto aInt = *(int64_t*)&A;
        auto bInt = *(int64_t*)&B;
#pragma GCC diagnostic pop

        // Different signs means they do not match.
        // Note: a negative floating point number is also negative as integer.
        if (std::signbit (aInt) != std::signbit (bInt))
            return false;

        // Find the difference in ULPs.
        return (std::abs (aInt - bInt) <= maxUlps);
    }

    inline bool almost_eq_zero(double A, double maxDiff = 1e-15) {
        return (std::abs (A) <= maxDiff);
    }
    FUNC_EQ(x, y);
    FUNC_EQ_ZERO(x);
    FUNC_LE(x, y);
    FUNC_LE_ZERO(x);
    FUNC_LT(x, y);
    FUNC_LT_ZERO(x);
    FUNC_GE(x, y);
    FUNC_GE_ZERO(x);
    FUNC_GT(x, y);
    FUNC_GT_ZERO(x);
}

namespace cmp = cmp_ulp;
/*

  Some
   _   _      _
  | | | | ___| |_ __   ___ _ __
  | |_| |/ _ \ | '_ \ / _ \ '__|
  |  _  |  __/ | |_) |  __/ |
  |_| |_|\___|_| .__/ \___|_|
             |_|

  functions
 */
namespace {
struct VectorLessX {
    bool operator() (Vector_t x1, Vector_t x2) {
        return cmp::lt (x1(0), x2(0));
    }
};

struct VectorLessY {
    bool operator() (Vector_t x1, Vector_t x2) {
        return cmp::lt(x1(1), x2 (1));
    }
};

struct VectorLessZ {
    bool operator() (Vector_t x1, Vector_t x2) {
        return cmp::lt(x1(2), x2(2));
    }
};

/**
   Calculate the maximum of coordinates of geometry,i.e the maximum of X,Y,Z
 */
Vector_t get_max_extent (std::vector<Vector_t>& coords) {
    const Vector_t x = *max_element (
        coords.begin (), coords.end (), VectorLessX ());
    const Vector_t y = *max_element (
        coords.begin (), coords.end (), VectorLessY ());
    const Vector_t z = *max_element (
        coords.begin (), coords.end (), VectorLessZ ());
    return Vector_t (x(0), y(1), z(2));
}


/*
   Compute the minimum of coordinates of geometry, i.e the minimum of X,Y,Z
 */
Vector_t get_min_extent (std::vector<Vector_t>& coords) {
    const Vector_t x = *min_element (
        coords.begin (), coords.end (), VectorLessX ());
    const Vector_t y = *min_element (
        coords.begin (), coords.end (), VectorLessY ());
    const Vector_t z = *min_element (
        coords.begin (), coords.end (), VectorLessZ ());
    return Vector_t (x(0), y(1), z(2));
}

/*
  write legacy VTK file of voxel mesh
*/
static void write_voxel_mesh (
    std::string fname,
    const std::unordered_map< int, std::unordered_set<int> >& ids,
    const Vector_t& hr_m,
    const Vektor<int,3>& nr,
    const Vector_t& origin
    ) {
    /*----------------------------------------------------------------------*/
    const size_t numpoints = 8 * ids.size ();
    std::ofstream of;

    *gmsg << level2 << "* Writing VTK file of voxel mesh '" << fname << "'" << endl;
    of.open (fname);
    PAssert (of.is_open ());
    of.precision (6);

    of << "# vtk DataFile Version 2.0" << std::endl;
    of << "generated using BoundaryGeometry::computeMeshVoxelization"
       << std::endl;
    of << "ASCII" << std::endl << std::endl;
    of << "DATASET UNSTRUCTURED_GRID" << std::endl;
    of << "POINTS " << numpoints << " float" << std::endl;

    const auto nr0_times_nr1 = nr[0] * nr[1];
    for (auto& elem: ids) {
        auto id = elem.first;
        int k = (id - 1) / nr0_times_nr1;
        int rest = (id - 1) % nr0_times_nr1;
        int j = rest / nr[0];
        int i = rest % nr[0];

        Vector_t P;
        P[0] = i * hr_m[0] + origin[0];
        P[1] = j * hr_m[1] + origin[1];
        P[2] = k * hr_m[2] + origin[2];

        of << P[0]           << " " << P[1]           << " " << P[2]           << std::endl;
        of << P[0] + hr_m[0] << " " << P[1]           << " " << P[2]           << std::endl;
        of << P[0]           << " " << P[1] + hr_m[1] << " " << P[2]           << std::endl;
        of << P[0] + hr_m[0] << " " << P[1] + hr_m[1] << " " << P[2]           << std::endl;
        of << P[0]           << " " << P[1]           << " " << P[2] + hr_m[2] << std::endl;
        of << P[0] + hr_m[0] << " " << P[1]           << " " << P[2] + hr_m[2] << std::endl;
        of << P[0]           << " " << P[1] + hr_m[1] << " " << P[2] + hr_m[2] << std::endl;
        of << P[0] + hr_m[0] << " " << P[1] + hr_m[1] << " " << P[2] + hr_m[2] << std::endl;
    }
    of << std::endl;
    const auto num_cells = ids.size ();
    of << "CELLS " << num_cells << " " << 9 * num_cells << std::endl;
    for (size_t i = 0; i < num_cells; i++)
        of << "8 "
           << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3 << " "
           << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << std::endl;
    of << "CELL_TYPES " << num_cells << std::endl;
    for (size_t i = 0; i <  num_cells; i++)
        of << "11" << std::endl;
    of << "CELL_DATA " << num_cells << std::endl;
    of << "SCALARS " << "cell_attribute_data" << " float " << "1" << std::endl;
    of << "LOOKUP_TABLE " << "default" << std::endl;
    for (size_t i = 0; i <  num_cells; i++)
        of << (float)(i) << std::endl;
    of << std::endl;
    of << "COLOR_SCALARS " << "BBoxColor " << 4 << std::endl;
    for (size_t i = 0; i < num_cells; i++) {
        of << "1.0" << " 1.0 " << "0.0 " << "1.0" << std::endl;
    }
    of << std::endl;
}
}

/*___________________________________________________________________________

  Triangle-cube intersection test.

  See:
  http://tog.acm.org/resources/GraphicsGems/gemsiii/triangleCube.c

 */


#define INSIDE 0
#define OUTSIDE 1

class Triangle {
public:
    Triangle () { }
    Triangle (const Vector_t& v1, const Vector_t& v2, const Vector_t& v3) {
        pts[0] = v1;
        pts[1] = v2;
        pts[2] = v3;
    }

    inline const Vector_t& v1() const {
        return pts[0];
    }
    inline double v1(int i) const {
        return pts[0][i];
    }
    inline const Vector_t& v2() const {
        return pts[1];
    }
    inline double v2(int i) const {
        return pts[1][i];
    }
    inline const Vector_t& v3() const {
        return pts[2];
    }
    inline double v3(int i) const {
        return pts[2][i];
    }


    inline void scale (
        const Vector_t& scaleby,
        const Vector_t& shiftby
        ) {
        pts[0][0] *= scaleby[0];
        pts[0][1] *= scaleby[1];
        pts[0][2] *= scaleby[2];
        pts[1][0] *= scaleby[0];
        pts[1][1] *= scaleby[1];
        pts[1][2] *= scaleby[2];
        pts[2][0] *= scaleby[0];
        pts[2][1] *= scaleby[1];
        pts[2][2] *= scaleby[2];
        pts[0] -= shiftby;
        pts[1] -= shiftby;
        pts[2] -= shiftby;
    }


    Vector_t pts[3];
};

/*___________________________________________________________________________*/

/* Which of the six face-plane(s) is point P outside of? */

static inline int
face_plane (
    const Vector_t& p
    ) {
    int outcode_fcmp = 0;

    if (cmp::gt(p[0],  0.5)) outcode_fcmp |= 0x01;
    if (cmp::lt(p[0], -0.5)) outcode_fcmp |= 0x02;
    if (cmp::gt(p[1],  0.5)) outcode_fcmp |= 0x04;
    if (cmp::lt(p[1], -0.5)) outcode_fcmp |= 0x08;
    if (cmp::gt(p[2],  0.5)) outcode_fcmp |= 0x10;
    if (cmp::lt(p[2], -0.5)) outcode_fcmp |= 0x20;

    return(outcode_fcmp);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . */

/* Which of the twelve edge plane(s) is point P outside of? */

static inline int
bevel_2d (
    const Vector_t& p
    ) {
    int outcode_fcmp = 0;

    if (cmp::gt( p[0] + p[1], 1.0)) outcode_fcmp |= 0x001;
    if (cmp::gt( p[0] - p[1], 1.0)) outcode_fcmp |= 0x002;
    if (cmp::gt(-p[0] + p[1], 1.0)) outcode_fcmp |= 0x004;
    if (cmp::gt(-p[0] - p[1], 1.0)) outcode_fcmp |= 0x008;
    if (cmp::gt( p[0] + p[2], 1.0)) outcode_fcmp |= 0x010;
    if (cmp::gt( p[0] - p[2], 1.0)) outcode_fcmp |= 0x020;
    if (cmp::gt(-p[0] + p[2], 1.0)) outcode_fcmp |= 0x040;
    if (cmp::gt(-p[0] - p[2], 1.0)) outcode_fcmp |= 0x080;
    if (cmp::gt( p[1] + p[2], 1.0)) outcode_fcmp |= 0x100;
    if (cmp::gt( p[1] - p[2], 1.0)) outcode_fcmp |= 0x200;
    if (cmp::gt(-p[1] + p[2], 1.0)) outcode_fcmp |= 0x400;
    if (cmp::gt(-p[1] - p[2], 1.0)) outcode_fcmp |= 0x800;

    return(outcode_fcmp);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  Which of the eight corner plane(s) is point P outside of?
*/
static inline int
bevel_3d (
    const Vector_t& p
    ) {
    int outcode_fcmp = 0;

    if (cmp::gt( p[0] + p[1] + p[2], 1.5)) outcode_fcmp |= 0x01;
    if (cmp::gt( p[0] + p[1] - p[2], 1.5)) outcode_fcmp |= 0x02;
    if (cmp::gt( p[0] - p[1] + p[2], 1.5)) outcode_fcmp |= 0x04;
    if (cmp::gt( p[0] - p[1] - p[2], 1.5)) outcode_fcmp |= 0x08;
    if (cmp::gt(-p[0] + p[1] + p[2], 1.5)) outcode_fcmp |= 0x10;
    if (cmp::gt(-p[0] + p[1] - p[2], 1.5)) outcode_fcmp |= 0x20;
    if (cmp::gt(-p[0] - p[1] + p[2], 1.5)) outcode_fcmp |= 0x40;
    if (cmp::gt(-p[0] - p[1] - p[2], 1.5)) outcode_fcmp |= 0x80;

    return(outcode_fcmp);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  Test the point "alpha" of the way from P1 to P2
  See if it is on a face of the cube
  Consider only faces in "mask"
*/

static inline int
check_point (
    const Vector_t& p1,
    const Vector_t& p2,
    const double alpha,
    const int mask
    ) {
    Vector_t plane_point;

#define LERP(a, b, t) (a + t * (b - a))
    // with C++20 we can use: std::lerp(a, b, t)
    plane_point[0] = LERP(p1[0], p2[0], alpha);
    plane_point[1] = LERP(p1[1], p2[1], alpha);
    plane_point[2] = LERP(p1[2], p2[2], alpha);
#undef LERP
    return(face_plane(plane_point) & mask);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  Compute intersection of P1 --> P2 line segment with face planes
  Then test intersection point to see if it is on cube face
  Consider only face planes in "outcode_diff"
  Note: Zero bits in "outcode_diff" means face line is outside of
*/
static inline int
check_line (
    const Vector_t& p1,
    const Vector_t& p2,
    const int outcode_diff
    ) {
    if ((0x01 & outcode_diff) != 0)
        if (check_point(p1,p2,( .5-p1[0])/(p2[0]-p1[0]),0x3e) == INSIDE) return(INSIDE);
    if ((0x02 & outcode_diff) != 0)
        if (check_point(p1,p2,(-.5-p1[0])/(p2[0]-p1[0]),0x3d) == INSIDE) return(INSIDE);
    if ((0x04 & outcode_diff) != 0)
        if (check_point(p1,p2,( .5-p1[1])/(p2[1]-p1[1]),0x3b) == INSIDE) return(INSIDE);
    if ((0x08 & outcode_diff) != 0)
        if (check_point(p1,p2,(-.5-p1[1])/(p2[1]-p1[1]),0x37) == INSIDE) return(INSIDE);
    if ((0x10 & outcode_diff) != 0)
        if (check_point(p1,p2,( .5-p1[2])/(p2[2]-p1[2]),0x2f) == INSIDE) return(INSIDE);
    if ((0x20 & outcode_diff) != 0)
        if (check_point(p1,p2,(-.5-p1[2])/(p2[2]-p1[2]),0x1f) == INSIDE) return(INSIDE);
    return(OUTSIDE);
}

/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  Test if 3D point is inside 3D triangle
*/
constexpr double EPS = 10e-15;
static inline int
SIGN3 (
    Vector_t A
    ) {
    return (((A[0] < EPS) ? 4 : 0) | ((A[0] > -EPS) ? 32 : 0) |
            ((A[1] < EPS) ? 2 : 0) | ((A[1] > -EPS) ? 16 : 0) |
            ((A[2] < EPS) ? 1 : 0) | ((A[2] > -EPS) ? 8  : 0));
}

static int
point_triangle_intersection (
    const Vector_t& p,
    const Triangle& t
    ) {
    /*
      First, a quick bounding-box test:
      If P is outside triangle bbox, there cannot be an intersection.
    */
    if (cmp::gt(p[0], std::max({t.v1(0), t.v2(0), t.v3(0)}))) return(OUTSIDE);
    if (cmp::gt(p[1], std::max({t.v1(1), t.v2(1), t.v3(1)}))) return(OUTSIDE);
    if (cmp::gt(p[2], std::max({t.v1(2), t.v2(2), t.v3(2)}))) return(OUTSIDE);
    if (cmp::lt(p[0], std::min({t.v1(0), t.v2(0), t.v3(0)}))) return(OUTSIDE);
    if (cmp::lt(p[1], std::min({t.v1(1), t.v2(1), t.v3(1)}))) return(OUTSIDE);
    if (cmp::lt(p[2], std::min({t.v1(2), t.v2(2), t.v3(2)}))) return(OUTSIDE);

    /*
      For each triangle side, make a vector out of it by subtracting vertexes;
      make another vector from one vertex to point P.
      The crossproduct of these two vectors is orthogonal to both and the
      signs of its X,Y,Z components indicate whether P was to the inside or
      to the outside of this triangle side.
    */
    const Vector_t vect12 = t.v1() - t.v2();
    const Vector_t vect1h = t.v1() - p;
    const Vector_t cross12_1p = cross (vect12, vect1h);
    const int sign12 = SIGN3(cross12_1p);      /* Extract X,Y,Z signs as 0..7 or 0...63 integer */

    const Vector_t vect23 = t.v2() - t.v3();
    const Vector_t vect2h = t.v2() - p;
    const Vector_t cross23_2p = cross (vect23, vect2h);
    const int sign23 = SIGN3(cross23_2p);

    const Vector_t vect31 = t.v3() - t.v1();
    const Vector_t vect3h = t.v3() - p;
    const Vector_t cross31_3p = cross (vect31, vect3h);
    const int sign31 = SIGN3(cross31_3p);

    /*
      If all three crossproduct vectors agree in their component signs,
      then the point must be inside all three.
      P cannot be OUTSIDE all three sides simultaneously.
    */
    return ((sign12 & sign23 & sign31) == 0) ? OUTSIDE : INSIDE;
}


/*. . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  This is the main algorithm procedure.
  Triangle t is compared with a unit cube,
  centered on the origin.
  It returns INSIDE (0) or OUTSIDE(1) if t
  intersects or does not intersect the cube.
*/
static int
triangle_intersects_cube (
    const Triangle& t
    ) {
    int v1_test;
    int v2_test;
    int v3_test;

    /*
      First compare all three vertexes with all six face-planes
      If any vertex is inside the cube, return immediately!
    */
    if ((v1_test = face_plane(t.v1())) == INSIDE) return(INSIDE);
    if ((v2_test = face_plane(t.v2())) == INSIDE) return(INSIDE);
    if ((v3_test = face_plane(t.v3())) == INSIDE) return(INSIDE);

   /*
     If all three vertexes were outside of one or more face-planes,
     return immediately with a trivial rejection!
   */
   if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);

   /*
     Now do the same trivial rejection test for the 12 edge planes
   */
   v1_test |= bevel_2d(t.v1()) << 8;
   v2_test |= bevel_2d(t.v2()) << 8;
   v3_test |= bevel_2d(t.v3()) << 8;
   if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);

   /*
     Now do the same trivial rejection test for the 8 corner planes
   */
   v1_test |= bevel_3d(t.v1()) << 24;
   v2_test |= bevel_3d(t.v2()) << 24;
   v3_test |= bevel_3d(t.v3()) << 24;
   if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);

   /*
     If vertex 1 and 2, as a pair, cannot be trivially rejected
     by the above tests, then see if the v1-->v2 triangle edge
     intersects the cube.  Do the same for v1-->v3 and v2-->v3./
     Pass to the intersection algorithm the "OR" of the outcode
     bits, so that only those cube faces which are spanned by
     each triangle edge need be tested.
   */
   if ((v1_test & v2_test) == 0)
       if (check_line (t.v1(), t.v2(), v1_test|v2_test) == INSIDE) return(INSIDE);
   if ((v1_test & v3_test) == 0)
       if (check_line (t.v1(), t.v3(), v1_test|v3_test) == INSIDE) return(INSIDE);
   if ((v2_test & v3_test) == 0)
       if (check_line (t.v2(), t.v3(), v2_test|v3_test) == INSIDE) return(INSIDE);

   /*
     By now, we know that the triangle is not off to any side,
     and that its sides do not penetrate the cube.  We must now
     test for the cube intersecting the interior of the triangle.
     We do this by looking for intersections between the cube
     diagonals and the triangle...first finding the intersection
     of the four diagonals with the plane of the triangle, and
     then if that intersection is inside the cube, pursuing
     whether the intersection point is inside the triangle itself.

     To find plane of the triangle, first perform crossproduct on
     two triangle side vectors to compute the normal vector.
   */
   Vector_t vect12 = t.v1() - t.v2();
   Vector_t vect13 = t.v1() - t.v3();
   Vector_t norm = cross (vect12, vect13);

   /*
     The normal vector "norm" X,Y,Z components are the coefficients
     of the triangles AX + BY + CZ + D = 0 plane equation.  If we
     solve the plane equation for X=Y=Z (a diagonal), we get
     -D/(A+B+C) as a metric of the distance from cube center to the
     diagonal/plane intersection.  If this is between -0.5 and 0.5,
     the intersection is inside the cube.  If so, we continue by
     doing a point/triangle intersection.
     Do this for all four diagonals.
   */
   double d = norm[0] * t.v1(0) + norm[1] * t.v1(1) + norm[2] * t.v1(2);

   /*
     if one of the diagonals is parallel to the plane, the other will
     intersect the plane
   */
   double denom = norm[0] + norm[1] + norm[2];
   if (cmp::eq_zero(std::abs(denom)) == false) {
       /* skip parallel diagonals to the plane; division by 0 can occure */
       Vector_t hitpp = d / denom;
       if (cmp::le(std::abs(hitpp[0]), 0.5))
           if (point_triangle_intersection(hitpp,t) == INSIDE)
               return(INSIDE);
   }
   denom = norm[0] + norm[1] - norm[2];
   if (cmp::eq_zero(std::abs(denom)) == false) {
       Vector_t hitpn;
       hitpn[2] = -(hitpn[0] = hitpn[1] = d / denom);
       if (cmp::le(std::abs(hitpn[0]), 0.5))
           if (point_triangle_intersection(hitpn,t) == INSIDE)
               return(INSIDE);
   }
   denom = norm[0] - norm[1] + norm[2];
   if (cmp::eq_zero(std::abs(denom)) == false) {
       Vector_t hitnp;
       hitnp[1] = -(hitnp[0] = hitnp[2] = d / denom);
       if (cmp::le(std::abs(hitnp[0]), 0.5))
           if (point_triangle_intersection(hitnp,t) == INSIDE)
               return(INSIDE);
   }
   denom = norm[0] - norm[1] - norm[2];
   if (cmp::eq_zero(std::abs(denom)) == false) {
       Vector_t hitnn;
       hitnn[1] = hitnn[2] = -(hitnn[0] = d / denom);
       if (cmp::le(std::abs(hitnn[0]), 0.5))
           if (point_triangle_intersection(hitnn,t) == INSIDE)
               return(INSIDE);
   }

   /*
     No edge touched the cube; no cube diagonal touched the triangle.
     We're done...there was no intersection.
   */
   return(OUTSIDE);
}

/*
 * Ray class, for use with the optimized ray-box intersection test
 * described in:
 *
 *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
 *      "An Efficient and Robust Ray-Box Intersection Algorithm"
 *      Journal of graphics tools, 10(1):49-54, 2005
 *
 */

class Ray {
public:
    Ray () { }
    Ray (Vector_t o, Vector_t d) {
        origin = o;
        direction = d;
        inv_direction = Vector_t (1/d[0], 1/d[1], 1/d[2]);
        sign[0] = (inv_direction[0] < 0);
        sign[1] = (inv_direction[1] < 0);
        sign[2] = (inv_direction[2] < 0);
    }
    Ray(const Ray &r) {
        origin = r.origin;
        direction = r.direction;
        inv_direction = r.inv_direction;
        sign[0] = r.sign[0]; sign[1] = r.sign[1]; sign[2] = r.sign[2];
    }
    const Ray &operator=(const Ray& a) = delete;

    Vector_t origin;
    Vector_t direction;
    Vector_t inv_direction;
    int sign[3];
};


/*
 * Axis-aligned bounding box class, for use with the optimized ray-box
 * intersection test described in:
 *
 *      Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
 *      "An Efficient and Robust Ray-Box Intersection Algorithm"
 *      Journal of graphics tools, 10(1):49-54, 2005
 *
 */

class Voxel {
public:
    Voxel () { }
    Voxel (const Vector_t &min, const Vector_t &max) {
        pts[0] = min;
        pts[1] = max;
    }
    inline void scale (
        const Vector_t& scale
        ) {
        pts[0][0] *= scale[0];
        pts[0][1] *= scale[1];
        pts[0][2] *= scale[2];
        pts[1][0] *= scale[0];
        pts[1][1] *= scale[1];
        pts[1][2] *= scale[2];
    }

    // (t0, t1) is the interval for valid hits
    bool intersect (
        const Ray& r,
        double& tmin,       // tmin and tmax are unchanged, if there is
        double& tmax        // no intersection
        ) const {
        double tmin_ = (pts[r.sign[0]][0]   - r.origin[0]) * r.inv_direction[0];
        double tmax_ = (pts[1-r.sign[0]][0] - r.origin[0]) * r.inv_direction[0];
        const double tymin = (pts[r.sign[1]][1]   - r.origin[1]) * r.inv_direction[1];
        const double tymax = (pts[1-r.sign[1]][1] - r.origin[1]) * r.inv_direction[1];
        if ( cmp::gt(tmin_, tymax) || cmp::gt(tymin, tmax_) )
            return 0;       // no intersection
        if (cmp::gt(tymin, tmin_))
            tmin_ = tymin;
        if (cmp::lt(tymax, tmax_))
            tmax_ = tymax;
        const double tzmin = (pts[r.sign[2]][2]   - r.origin[2]) * r.inv_direction[2];
        const double tzmax = (pts[1-r.sign[2]][2] - r.origin[2]) * r.inv_direction[2];
        if ( cmp::gt(tmin_, tzmax) || cmp::gt(tzmin, tmax_) )
            return 0;       // no intersection
        if (cmp::gt(tzmin, tmin_))
            tmin_ = tzmin;
        tmin = tmin_;
        if (cmp::lt(tzmax, tmax_))
            tmax_ = tzmax;
        tmax = tmax_;
        return cmp::ge_zero(tmax);
    }

    inline bool intersect (
        const Ray& r
        ) const {
        double tmin = 0.0;
        double tmax = 0.0;
        return intersect(r, tmin, tmax);
    }

    inline int intersect (
        const Triangle& t
        ) const {
        Voxel v_ = *this;
        Triangle t_ = t;
        const Vector_t scaleby = 1.0 / v_.extent();
        v_.scale (scaleby);
        t_.scale (scaleby , v_.pts[0] + 0.5);
        return triangle_intersects_cube (t_);
    }

    inline Vector_t extent () const {
        return (pts[1] - pts[0]);
    }

    inline bool isInside  (
        const Vector_t& P
        ) const {
        return (
                cmp::ge(P[0], pts[0][0])
                && cmp::ge(P[1], pts[0][1])
                && cmp::ge(P[2], pts[0][2])
                && cmp::le(P[0], pts[1][0])
                && cmp::le(P[1], pts[1][1])
                && cmp::le(P[2], pts[1][2]));
    }

    Vector_t pts[2];
};

static inline Vector_t normalVector (
    const Vector_t& A,
    const Vector_t& B,
    const Vector_t& C
    ) {
    const Vector_t N = cross (B - A, C - A);
    const double magnitude = std::sqrt (SQR (N (0)) + SQR (N (1)) + SQR (N (2)));
    PAssert (cmp::gt_zero(magnitude)); // in case we have degenerated triangles
    return N / magnitude;
}

// Calculate the area of triangle given by id.
static inline double computeArea (
    const Vector_t& A,
    const Vector_t& B,
    const Vector_t& C
    ) {
    const Vector_t AB = A - B;
    const Vector_t AC = C - A;
    return(0.5 * std::sqrt (dot (AB, AB) * dot (AC, AC) - dot (AB, AC) * dot (AB, AC)));
}


/*
  ____                           _
 / ___| ___  ___  _ __ ___   ___| |_ _ __ _   _
| |  _ / _ \/ _ \| '_ ` _ \ / _ \ __| '__| | | |
| |_| |  __/ (_) | | | | | |  __/ |_| |  | |_| |
 \____|\___|\___/|_| |_| |_|\___|\__|_|   \__, |
                                          |___/
*/

BoundaryGeometry::BoundaryGeometry() :
    Definition (
        SIZE, "GEOMETRY", "The \"GEOMETRY\" statement defines the beam pipe geometry.") {

    itsAttr[FGEOM] = Attributes::makeString
        ("FGEOM",
         "Specifies the geometry file [H5hut]",
         "");

    itsAttr[TOPO] = Attributes::makePredefinedString
        ("TOPO",
         "If FGEOM is selected topo is over-written. ",
         {"RECTANGULAR", "BOXCORNER", "ELLIPTIC"},
         "ELLIPTIC");

    itsAttr[LENGTH] = Attributes::makeReal
        ("LENGTH",
         "Specifies the length of a tube shaped elliptic beam pipe [m]",
         1.0);

    itsAttr[S] = Attributes::makeReal
        ("S",
         "Specifies the start of a tube shaped elliptic beam pipe [m]",
         0.0);

    itsAttr[A] = Attributes::makeReal
        ("A",
         "Specifies the major semi-axis of a tube shaped elliptic beam pipe [m]",
         0.025);

    itsAttr[B] = Attributes::makeReal
        ("B",
         "Specifies the major semi-axis of a tube shaped elliptic beam pipe [m]",
         0.025);

    itsAttr[L1] = Attributes::makeReal
        ("L1",
         "In case of BOXCORNER Specifies first part with height == B [m]",
         0.5);

    itsAttr[L2] = Attributes::makeReal
        ("L2",
         "In case of BOXCORNER Specifies first second with height == B-C [m]",
         0.2);

    itsAttr[C] = Attributes::makeReal
        ("C",
         "In case of BOXCORNER Specifies height of corner C [m]",
         0.01);

    itsAttr[XYZSCALE] = Attributes::makeReal
        ("XYZSCALE",
         "Multiplicative scaling factor for coordinates ",
         1.0);

    itsAttr[XSCALE] = Attributes::makeReal
        ("XSCALE",
         "Multiplicative scaling factor for X coordinates ",
         1.0);

    itsAttr[YSCALE] = Attributes::makeReal
        ("YSCALE",
         "Multiplicative scaling factor for Y coordinates ",
         1.0);

    itsAttr[ZSCALE] = Attributes::makeReal
        ("ZSCALE",
         "Multiplicative scaling factor for Z coordinates ",
         1.0);

    itsAttr[ZSHIFT] = Attributes::makeReal
        ("ZSHIFT",
         "Shift in z direction",
         0.0);

    itsAttr[INSIDEPOINT] = Attributes::makeRealArray
        ("INSIDEPOINT", "A point inside the geometry");

    registerOwnership(AttributeHandler::STATEMENT);

    BoundaryGeometry* defGeometry = clone ("UNNAMED_GEOMETRY");
    defGeometry->builtin = true;

    Tinitialize_m =   IpplTimings::getTimer ("Initialize geometry");
    TisInside_m =     IpplTimings::getTimer ("Inside test");
    TfastIsInside_m = IpplTimings::getTimer ("Fast inside test");
    TRayTrace_m =     IpplTimings::getTimer ("Ray tracing");
    TPartInside_m =   IpplTimings::getTimer ("Particle Inside");

    h5FileName_m = Attributes::getString (itsAttr[FGEOM]);

    try {
        defGeometry->update ();
        OpalData::getInstance ()->define (defGeometry);
    } catch (...) {
        delete defGeometry;
    }
    gsl_rng_env_setup();
    randGen_m = gsl_rng_alloc(gsl_rng_default);

    if (!h5FileName_m.empty ())
        initialize ();
}

BoundaryGeometry::BoundaryGeometry(
    const std::string& name,
    BoundaryGeometry* parent
    ) : Definition (name, parent) {
    gsl_rng_env_setup();
    randGen_m = gsl_rng_alloc(gsl_rng_default);

    Tinitialize_m =   IpplTimings::getTimer ("Initialize geometry");
    TisInside_m =     IpplTimings::getTimer ("Inside test");
    TfastIsInside_m = IpplTimings::getTimer ("Fast inside test");
    TRayTrace_m =     IpplTimings::getTimer ("Ray tracing");
    TPartInside_m =   IpplTimings::getTimer ("Particle Inside");

    h5FileName_m = Attributes::getString (itsAttr[FGEOM]);
    if (!h5FileName_m.empty ())
        initialize ();
 }

BoundaryGeometry::~BoundaryGeometry() {
    gsl_rng_free(randGen_m);
}

bool BoundaryGeometry::canReplaceBy (Object* object) {
    // Can replace only by another GEOMETRY.
    return dynamic_cast<BGeometryBase*>(object) != 0;
}

BoundaryGeometry* BoundaryGeometry::clone (const std::string& name) {
    return new BoundaryGeometry (name, this);
}

void BoundaryGeometry::update () {
    if (getOpalName ().empty ()) setOpalName ("UNNAMED_GEOMETRY");
}


void BoundaryGeometry::execute () {
    update ();
    Tinitialize_m =   IpplTimings::getTimer ("Initialize geometry");
    TisInside_m =     IpplTimings::getTimer ("Inside test");
    TfastIsInside_m = IpplTimings::getTimer ("Fast inside test");
    TRayTrace_m =     IpplTimings::getTimer ("Ray tracing");
    TPartInside_m =   IpplTimings::getTimer ("Particle Inside");
}

BoundaryGeometry* BoundaryGeometry::find (const std::string& name) {
    BoundaryGeometry* geom = dynamic_cast<BoundaryGeometry*>(
        OpalData::getInstance ()->find (name));

    if (geom == 0)
        throw OpalException ("BoundaryGeometry::find()", "Geometry \""
                             + name + "\" not found.");
    return geom;
}

void BoundaryGeometry::updateElement (ElementBase* /*element*/) {
}

int
BoundaryGeometry::intersectTriangleVoxel (
    const int triangle_id,
    const int i,
    const int j,
    const int k
    ) {
    const Triangle t(
        getPoint (triangle_id, 1),
        getPoint (triangle_id, 2),
        getPoint (triangle_id, 3)
        );

    const Vector_t P(
        i * voxelMesh_m.sizeOfVoxel [0] + voxelMesh_m.minExtent[0],
        j * voxelMesh_m.sizeOfVoxel [1] + voxelMesh_m.minExtent[1],
        k * voxelMesh_m.sizeOfVoxel [2] + voxelMesh_m.minExtent[2]
        );

    Voxel v(P, P+voxelMesh_m.sizeOfVoxel);

    return v.intersect (t);
}

/*
  Find the 3D intersection of a line segment, ray or line with a triangle.

  Input:
        kind: type of test: SEGMENT, RAY or LINE
        P0, P0: defining
            a line segment from P0 to P1 or
            a ray starting at P0 with directional vector P1-P0 or
            a line through P0 and P1
        V0, V1, V2: the triangle vertices

  Output:
        I: intersection point (when it exists)

  Return values for line segment and ray test :
        -1 = triangle is degenerated (a segment or point)
        0 =  disjoint (no intersect)
        1 =  are in the same plane
        2 =  intersect in unique point I1

  Return values for line intersection test :
        -1: triangle is degenerated (a segment or point)
        0:  disjoint (no intersect)
        1:  are in the same plane
        2:  intersect in unique point I1, with r < 0.0
        3:  intersect in unique point I1, with 0.0 <= r <= 1.0
        4:  intersect in unique point I1, with 1.0 < r

  For algorithm and implementation see:
  http://geomalgorithms.com/a06-_intersect-2.html

  Copyright 2001 softSurfer, 2012 Dan Sunday
  This code may be freely used and modified for any purpose
  providing that this copyright notice is included with it.
  SoftSurfer makes no warranty for this code, and cannot be held
  liable for any real or imagined damage resulting from its use.
  Users of this code must verify correctness for their application.
 */

int
BoundaryGeometry::intersectLineTriangle (
    const enum INTERSECTION_TESTS kind,
    const Vector_t& P0,
    const Vector_t& P1,
    const int triangle_id,
    Vector_t& I
    ) {
    const Vector_t V0 = getPoint (triangle_id, 1);
    const Vector_t V1 = getPoint (triangle_id, 2);
    const Vector_t V2 = getPoint (triangle_id, 3);

    // get triangle edge vectors and plane normal
    const Vector_t u = V1 - V0;         // triangle vectors
    const Vector_t v = V2 - V0;
    const Vector_t n = cross (u, v);
    if (n == (Vector_t)0)               // triangle is degenerate
        return -1;                      // do not deal with this case

    const Vector_t dir = P1 - P0;       // ray direction vector
    const Vector_t w0 = P0 - V0;
    const double a = -dot(n,w0);
    const double b = dot(n,dir);
    if (cmp::eq_zero(b)) {              // ray is  parallel to triangle plane
        if (cmp::eq_zero(a)) {          // ray lies in triangle plane
            return 1;
        } else {                        // ray disjoint from plane
            return 0;
        }
    }

    // get intersect point of ray with triangle plane
    const double r = a / b;
    switch (kind) {
    case RAY:
        if (cmp::lt_zero(r)) {          // ray goes away from triangle
            return 0;                   // => no intersect
        }
        break;
    case SEGMENT:
        if (cmp::lt_zero(r) || cmp::lt(1.0, r)) { // intersection on extended
            return 0;                   // segment
        }
        break;
    case LINE:
        break;
    };
    I = P0 + r * dir;                   // intersect point of ray and plane

    // is I inside T?
    const double uu = dot(u,u);
    const double uv = dot(u,v);
    const double vv = dot(v,v);
    const Vector_t w = I - V0;
    const double wu = dot(w,u);
    const double wv = dot(w,v);
    const double D = uv * uv - uu * vv;

    // get and test parametric coords
    const double s = (uv * wv - vv * wu) / D;
    if (cmp::lt_zero(s) || cmp::gt(s, 1.0)) {           // I is outside T
        return 0;
    }
    const double t = (uv * wu - uu * wv) / D;
    if (cmp::lt_zero(t) || cmp::gt((s + t), 1.0)) {     // I is outside T
        return 0;
    }
    // intersection point is in triangle
    if (cmp::lt_zero(r)) {              // in extended segment in opposite
        return 2;                       // direction of ray
    } else if (cmp::ge_zero(r) && cmp::le(r, 1.0)) { // in segment
        return 3;
    } else {                            // in extended segment in
        return 4;                       // direction of ray
    }
}

static inline double magnitude (
    const Vector_t& v
    ) {
    return std::sqrt (dot (v,v));
}

bool
BoundaryGeometry::isInside (
    const Vector_t& P                    // [in] pt to test
    ) {

    /*
      select a "close" reference pt outside the bounding box
    */
    // right boundary of bounding box (x direction)
    double x = minExtent_m[0] - 0.01;
    double distance = P[0] - x;
    Vector_t ref_pt {x, P[1], P[2]};

    // left boundary of bounding box (x direction)
    x = maxExtent_m[0] + 0.01;
    if (cmp::lt(x - P[0], distance)) {
        distance = x - P[0];
        ref_pt = {x, P[1], P[2]};
    }

    // lower boundary of bounding box (y direction)
    double y = minExtent_m[1] - 0.01;
    if (cmp::lt(P[1] - y, distance)) {
        distance = P[1] -y;
        ref_pt = {P[0], y, P[1]};
    }

    // upper boundary of bounding box (y direction)
    y = maxExtent_m[1] + 0.01;
    if (cmp::lt(y - P[1], distance)) {
        distance = y - P[1];
        ref_pt = {P[0], y, P[2]};
    }
    // front boundary of bounding box (z direction)
    double z = minExtent_m[2] - 0.01;
    if (cmp::lt(P[2] - z, distance)) {
        distance = P[2] - z;
        ref_pt = {P[0], P[1], z};
    }
    // back boundary of bounding box (z direction)
    z = maxExtent_m[2] + 0.01;
    if (cmp::lt(z - P[2], distance)) {
        ref_pt = {P[0], P[1], z};
    }

    /*
      the test returns the number of intersections =>
      since the reference point is outside, P is inside
      if the result is odd.
    */
    int k = fastIsInside (ref_pt, P);
    return (k % 2) == 1;
}

/*
  searching a point inside the geometry.

  sketch of the algorithm:
  In a first step, we try to find a line segment defined by one
  point outside the bounding box and a point somewhere inside the
  bounding box which has intersects with the geometry.

  If the number of intersections is odd, the center point is inside
  the geometry and we are already done.

  If the number of intersections is even, there must be points on
  this line segment which are inside the geometry. In the next step
  we have to find one if these points.


  A bit more in detail:

  1. Finding a line segment intersecting the geometry
  For the fast isInside test it is of advantage to choose line segments
  parallel to the X, Y or Z axis. In this implementation we choose as
  point outside the bounding box a point on an axis but close to the
  bounding box and the center of the bounding box. This gives us six
  line segments to test. This covers not all possible geometries but
  most likely almost all. If not, it's easy to extend.

  2. Searching for a point inside the geometry
  In the first step we get a line segment from which we know, that one
  point is ouside the geometry (P_out) and the other inside the bounding
  box (Q). We also know the number of intersections n_i of this line
  segment with the geometry.

  If n_i is odd, Q is inside the boundary!

  while (true); do
      bisect the line segment [P_out, Q], let B the bisecting point.

      compute number of intersections of the line segment [P_out, B]
      and the geometry.

      If the number of intersections is odd, then B is inside the geometry
      and we are done. Set P_in = B and exit loop.

      Otherwise we have either no or an even number of intersections.
      In both cases this implies that B is a point outside the geometry.

      If the number of intersection of [P_out, B] is even but not equal zero,
      it might be that *all* intersections are in this line segment and none in
      [B, Q].
      In this case we continue with the line segment [P_out, Q] = [P_out, B],
      otherwise with the line segment [P_out, Q] = [B, Q].
*/
bool
BoundaryGeometry::findInsidePoint (
    void
    ) {
    *gmsg << level2 << "* Searching for a point inside the geometry..." << endl;
    /*
      find line segment
    */
    Vector_t Q {(maxExtent_m + minExtent_m) / 2};
    std::vector<Vector_t> P_outs {
        {minExtent_m[0]-0.01, Q[1], Q[2]},
        {maxExtent_m[0]+0.01, Q[1], Q[2]},
        {Q[0], minExtent_m[1]-0.01, Q[2]},
        {Q[0], maxExtent_m[1]+0.01, Q[2]},
        {Q[0], Q[1], minExtent_m[2]-0.01},
        {Q[0], Q[1], maxExtent_m[2]+0.01}
    };
    int n_i = 0;
    Vector_t P_out;
    for (const auto& P: P_outs) {
        n_i = fastIsInside (P, Q);
        if (n_i != 0) {
            P_out = P;
            break;
        }
    }
    if (n_i == 0) {
        // this is possible with some obscure geometries.
        return false;
    }

    /*
      if the number of intersections is odd, Q is inside the geometry
    */
    if (n_i % 2 == 1) {
        insidePoint_m = Q;
        return true;
    }
    while (true) {
        Vector_t B {(P_out + Q) / 2};
        int n = fastIsInside (P_out, B);
        if (n % 2 == 1) {
            insidePoint_m = B;
            return true;
        } else if (n == n_i) {
            Q = B;
        } else {
            P_out = B;
        }
        n_i = n;
    }
    // never reached
    return false;
}

/*
  Game plan:
  Count number of intersection of the line segment defined by P and a reference
  pt with the boundary. If the reference pt is inside the boundary and the number
  of intersections is even, then P is inside the geometry. Otherwise P is outside.
  To count the number of intersection, we divide the line segment in N segments
  and run the line-segment boundary intersection test for all these segments.
  N must be choosen carefully. It shouldn't be to large to avoid needless test.
 */
int
BoundaryGeometry::fastIsInside (
    const Vector_t& reference_pt,        // [in] reference pt inside the boundary
    const Vector_t& P                    // [in] pt to test
    ) {
    const Voxel c (minExtent_m, maxExtent_m);
    if (!c.isInside (P)) return 1;
    IpplTimings::startTimer (TfastIsInside_m);
#ifdef ENABLE_DEBUG
    int saved_flags = debugFlags_m;
    if (debugFlags_m & debug_fastIsInside) {
        *gmsg << "* " << __func__ << ": "
              << "reference_pt=" << reference_pt
              << ",  P=" << P << endl;
        debugFlags_m |= debug_intersectTinyLineSegmentBoundary;
    }
#endif
    const Vector_t v = reference_pt - P;
    const int N = std::ceil (magnitude (v) / std::min ({voxelMesh_m.sizeOfVoxel [0],
                                              voxelMesh_m.sizeOfVoxel [1],
                                              voxelMesh_m.sizeOfVoxel [2]}));
    const Vector_t v_ = v / N;
    Vector_t P0 = P;
    Vector_t P1 = P + v_;
    Vector_t I;
    int triangle_id = -1;
    int result = 0;
    for (int i = 0; i < N; i++) {
        result += intersectTinyLineSegmentBoundary (P0, P1, I, triangle_id);
        P0 = P1;
        P1 += v_;
    }
#ifdef ENABLE_DEBUG
    if (debugFlags_m & debug_fastIsInside) {
        *gmsg << "* " << __func__ << ": "
              << "result: " << result << endl;
        debugFlags_m = saved_flags;
    }
#endif
    IpplTimings::stopTimer (TfastIsInside_m);
    return result;
}

/*
  P must be *inside* the boundary geometry!

  return value:
    0   no intersection
    1   intersection found, I is set to the first intersection coordinates in
        ray direction
 */
int
BoundaryGeometry::intersectRayBoundary (
    const Vector_t& P,
    const Vector_t& v,
    Vector_t& I
    ) {
    IpplTimings::startTimer (TRayTrace_m);
#ifdef ENABLE_DEBUG
    int saved_flags = debugFlags_m;
    if (debugFlags_m & debug_intersectRayBoundary) {
        *gmsg << "* " << __func__ << ": "
              << "  ray: "
              << "  origin=" << P
              << "  dir=" << v
              << endl;
        debugFlags_m |= debug_intersectLineSegmentBoundary;
    }
#endif
    /*
      set P1 to intersection of ray with bbox of voxel mesh
      run line segment boundary intersection test with P and P1
     */
    Ray r = Ray (P, v);
    Voxel c = Voxel (voxelMesh_m.minExtent+0.25*voxelMesh_m.sizeOfVoxel,
                     voxelMesh_m.maxExtent-0.25*voxelMesh_m.sizeOfVoxel);
    double tmin = 0.0;
    double tmax = 0.0;
    c.intersect (r, tmin, tmax);
    int triangle_id = -1;
    int result = (intersectLineSegmentBoundary (
                      P, P + (tmax*v),
                      I, triangle_id) > 0) ? 1 : 0;
#ifdef ENABLE_DEBUG
    if (debugFlags_m & debug_intersectRayBoundary) {
        *gmsg << "* " << __func__ << ": "
              << "  result=" << result
              << "  I=" << I
              << endl;
        debugFlags_m = saved_flags;
    }
#endif
    IpplTimings::stopTimer (TRayTrace_m);
    return result;
}

/*
  Map point to unique voxel ID.

  Remember:
  * hr_m:  is the  mesh size
  * nr_m:  number of mesh points
  */
inline int
BoundaryGeometry::mapVoxelIndices2ID (
    const int i,
    const int j,
    const int k
    ) {
    if (i < 0 || i >= voxelMesh_m.nr_m[0] ||
        j < 0 || j >= voxelMesh_m.nr_m[1] ||
        k < 0 || k >= voxelMesh_m.nr_m[2]) {
        return 0;
    }
    return 1 + k * voxelMesh_m.nr_m[0] * voxelMesh_m.nr_m[1] + j * voxelMesh_m.nr_m[0] + i;
}

#define mapPoint2VoxelIndices(pt, i, j, k) {                            \
        i = floor ((pt[0] - voxelMesh_m.minExtent [0]) / voxelMesh_m.sizeOfVoxel[0]); \
        j = floor ((pt[1] - voxelMesh_m.minExtent [1]) / voxelMesh_m.sizeOfVoxel[1]); \
        k = floor ((pt[2] - voxelMesh_m.minExtent [2]) / voxelMesh_m.sizeOfVoxel[2]); \
        if (!(0 <= i && i < voxelMesh_m.nr_m[0] &&                      \
              0 <= j && j < voxelMesh_m.nr_m[1] &&                      \
              0 <= k && k < voxelMesh_m.nr_m[2])) {                     \
            *gmsg << level2                                             \
                  << "* " << __func__ << ":"                            \
                  << "  WARNING: pt=" << pt                             \
                  << "  is outside the bbox"                            \
                  << "  i=" << i                                        \
                  << "  j=" << j                                        \
                  << "  k=" << k                                        \
                  << endl;                                              \
        }                                                               \
    }

inline Vector_t
BoundaryGeometry::mapIndices2Voxel (
    const int i,
    const int j,
    const int k
    ) {
    return Vector_t (
        i * voxelMesh_m.sizeOfVoxel [0] + voxelMesh_m.minExtent[0],
        j * voxelMesh_m.sizeOfVoxel [1] + voxelMesh_m.minExtent[1],
        k * voxelMesh_m.sizeOfVoxel [2] + voxelMesh_m.minExtent[2]);
}

inline Vector_t
BoundaryGeometry::mapPoint2Voxel (
    const Vector_t& pt
    ) {
    const int i = std::floor ((pt[0] - voxelMesh_m.minExtent [0]) / voxelMesh_m.sizeOfVoxel [0]);
    const int j = std::floor ((pt[1] - voxelMesh_m.minExtent [1]) / voxelMesh_m.sizeOfVoxel [1]);
    const int k = std::floor ((pt[2] - voxelMesh_m.minExtent [2]) / voxelMesh_m.sizeOfVoxel [2]);

    return mapIndices2Voxel (i, j, k);
}


inline void
BoundaryGeometry::computeMeshVoxelization (void) {

    for (unsigned int triangle_id = 0; triangle_id < Triangles_m.size(); triangle_id++) {
        Vector_t v1 = getPoint (triangle_id, 1);
        Vector_t v2 = getPoint (triangle_id, 2);
        Vector_t v3 = getPoint (triangle_id, 3);
        Vector_t bbox_min = {
            std::min({v1[0], v2[0], v3[0]}),
            std::min({v1[1], v2[1], v3[1]}),
            std::min({v1[2], v2[2], v3[2]}) };
        Vector_t bbox_max = {
            std::max({v1[0], v2[0], v3[0]}),
            std::max({v1[1], v2[1], v3[1]}),
            std::max({v1[2], v2[2], v3[2]}) };
        int i_min, j_min, k_min;
        int i_max, j_max, k_max;
        mapPoint2VoxelIndices (bbox_min, i_min, j_min, k_min);
        mapPoint2VoxelIndices (bbox_max, i_max, j_max, k_max);

        for (int i = i_min; i <= i_max; i++) {
            for (int j = j_min; j <= j_max; j++) {
                for (int k = k_min; k <= k_max; k++) {
                    // test if voxel (i,j,k) has an intersection with triangle
                    if (intersectTriangleVoxel (triangle_id, i, j, k) == INSIDE) {
                        auto id = mapVoxelIndices2ID (i, j, k);
                        voxelMesh_m.ids [id].insert (triangle_id);
                    }
                }
            }
        }
    } // for_each triangle
    *gmsg << level2 << "* Mesh voxelization done" << endl;

    // write voxel mesh into VTK file
    if (Ippl::myNode() == 0 && Options::enableVTK) {
        std::string vtkFileName = Util::combineFilePath({
            OpalData::getInstance()->getAuxiliaryOutputDirectory(),
            "testBBox.vtk"
        });
        bool writeVTK = false;

        if (!boost::filesystem::exists(vtkFileName)) {
            writeVTK = true;
        } else {
            std::time_t t_geom = boost::filesystem::last_write_time(h5FileName_m);
            std::time_t t_vtk = boost::filesystem::last_write_time(vtkFileName);
            if (std::difftime(t_geom,t_vtk) > 0)
                writeVTK = true;
        }

        if (writeVTK) {
            write_voxel_mesh (vtkFileName,
                              voxelMesh_m.ids,
                              voxelMesh_m.sizeOfVoxel,
                              voxelMesh_m.nr_m,
                              voxelMesh_m.minExtent);
        }
    }
}

void BoundaryGeometry::initialize () {

    class Local {

    public:

        static void computeGeometryInterval (BoundaryGeometry* bg) {

            bg->minExtent_m = get_min_extent (bg->Points_m);
            bg->maxExtent_m = get_max_extent (bg->Points_m);

            /*
              Calculate the maximum size of triangles. This value will be used to
              define the voxel size
            */
            double longest_edge_max_m = 0.0;
            for (unsigned int i = 0; i < bg->Triangles_m.size(); i++) {
                // compute length of longest edge
                const Vector_t x1 = bg->getPoint (i, 1);
                const Vector_t x2 = bg->getPoint (i, 2);
                const Vector_t x3 = bg->getPoint (i, 3);
                const double length_edge1 = std::sqrt (
                    SQR (x1[0] - x2[0]) + SQR (x1[1] - x2[1]) + SQR (x1[2] - x2[2]));
                const double length_edge2 = std::sqrt (
                    SQR (x3[0] - x2[0]) + SQR (x3[1] - x2[1]) + SQR (x3[2] - x2[2]));
                const double length_edge3 = std::sqrt (
                    SQR (x3[0] - x1[0]) + SQR (x3[1] - x1[1]) + SQR (x3[2] - x1[2]));

                double max = length_edge1;
                if (length_edge2 > max) max = length_edge2;
                if (length_edge3 > max) max = length_edge3;

                // save min and max of length of longest edge
                if (longest_edge_max_m < max) longest_edge_max_m = max;
            }

            /*
              In principal the number of discretization nr_m is the extent of
              the geometry divided by the extent of the largest triangle. Whereby
              the extent of a triangle is defined as the lenght of its longest
              edge. Thus the largest triangle is the triangle with the longest edge.

              But if the hot spot, i.e., the multipacting/field emission zone is
              too small that the normal bounding box covers the whole hot spot, the
              expensive triangle-line intersection tests will be frequently called.
              In these cases, we have to use smaller bounding box size to speed up
              simulation.

              Todo:
              The relation between bounding box size and simulation time step &
              geometry shape maybe need to be summarized and modeled in a more
              flexible manner and could be adjusted in input file.
            */
            Vector_t extent = bg->maxExtent_m - bg->minExtent_m;
            bg->voxelMesh_m.nr_m (0) = 16 * (int)std::floor (extent [0] / longest_edge_max_m);
            bg->voxelMesh_m.nr_m (1) = 16 * (int)std::floor (extent [1] / longest_edge_max_m);
            bg->voxelMesh_m.nr_m (2) = 16 * (int)std::floor (extent [2] / longest_edge_max_m);

            bg->voxelMesh_m.sizeOfVoxel = extent / bg->voxelMesh_m.nr_m;
            bg->voxelMesh_m.minExtent = bg->minExtent_m - 0.5 * bg->voxelMesh_m.sizeOfVoxel;
            bg->voxelMesh_m.maxExtent = bg->maxExtent_m + 0.5 * bg->voxelMesh_m.sizeOfVoxel;
            bg->voxelMesh_m.nr_m += 1;
        }

        /*
          To speed up ray-triangle intersection tests, the normal vector of
          all triangles are pointing inward. Since this is clearly not
          guaranteed for the triangles in the H5hut file, this must be checked
          for each triangle and - if necessary changed - after reading the mesh.

          To test whether the normal of a triangle is pointing inward or outward,
          we choose a random point P close to the center of the triangle and test
          whether this point is inside or outside the geometry. The way we choose
          P guarantees that the vector spanned by P and a vertex of the triangle
          points into the same direction as the normal vector. From this it
          follows that if P is inside the geometry the normal vector is pointing
          to the inside and vise versa.

          Since the inside-test is computational expensive we perform this test
          for one reference triangle T (per sub-mesh) only. Knowing the adjacent
          triangles for all three edges of a triangle for all triangles of the
          mesh facilitates another approach using the orientation of the
          reference triangle T. Assuming that the normal vector of T points to
          the inside of the geometry an adjacent triangle of T has an inward
          pointing normal vector if and only if it has the same orientation as
          T.

          Starting with the reference triangle T we can change the orientation
          of the adjancent triangle of T and so on.

          NOTE: For the time being we do not make use of the inward pointing
          normals.

          FIXME: Describe the basic ideas behind the following comment! Without
          it is completely unclear.

          Following combinations are possible:
              1,1 && 2,2   1,2 && 2,1   1,3 && 2,1
              1,1 && 2,3   1,2 && 2,3   1,3 && 2,2
              1,1 && 3,2   1,2 && 3,1   1,3 && 3,1
              1,1 && 3,3   1,2 && 3,3   1,3 && 3,2

             (2,1 && 1,2) (2,2 && 1,1) (2,3 && 1,1)
             (2,1 && 1,3) (2,2 && 1,3) (2,3 && 1,2)
              2,1 && 3,2   2,2 && 3,1   2,3 && 3,1
              2,1 && 3,3   2,2 && 3,3   2,3 && 3,2

             (3,1 && 1,2) (3,2 && 1,1) (3,3 && 1,1)
             (3,1 && 1,3) (3,2 && 1,3) (3,3 && 1,2)
             (3,1 && 2,2) (3,2 && 2,1) (3,3 && 2,1)
             (3,1 && 2,3) (3,2 && 2,3) (3,3 && 2,2)

          Note:
          Since we find vertices with lower enumeration first, we
          can ignore combinations in ()

                  2 2           2 3           3 2           3 3
                   *             *             *             *
                  /|\           /|\           /|\           /|\
                 / | \         / | \         / | \         / | \
                /  |  \       /  |  \       /  |  \       /  |  \
               /   |   \     /   |   \     /   |   \     /   |   \
              *----*----*   *----*----*   *----*----*   *----*----*
              3   1 1   3   3   1 1   2   2   1 1   3   2   1 1   2
diff:            (1,1)         (1,2)         (2,1)         (2,2)
change orient.:   yes           no            no            yes


                  2 1           2 3           3 1           3 3
                   *             *             *             *
                  /|\           /|\           /|\           /|\
                 / | \         / | \         / | \         / | \
                /  |  \       /  |  \       /  |  \       /  |  \
               /   |   \     /   |   \     /   |   \     /   |   \
              *----*----*   *----*----*   *----*----*   *----*----*
              3   1 2   3   3   1 2   1   2   1 2   3   2   1 2   1
diff:            (1,-1)        (1,1)         (2,-1)        (2,1)
change orient.:   no            yes           yes           no


                  2 1           2 2           3 1           3 2
                   *             *             *             *
                  /|\           /|\           /|\           /|\
                 / | \         / | \         / | \         / | \
                /  |  \       /  |  \       /  |  \       /  |  \
               /   |   \     /   |   \     /   |   \     /   |   \
              *----*----*   *----*----*   *----*----*   *----*----*
              3   1 3   2   3   1 3   1   2   1 3   2   2   1 3   1
diff:            (1,-2)        (1,-1)        (2,-2)        (2,-1)
change orient.:   yes           no            no            yes

                                              3 2           3 3
                                               *             *
                                              /|\           /|\
                                             / | \         / | \
                                            /  |  \       /  |  \
                                           /   |   \     /   |   \
                                          *----*----*   *----*----*
                                          1   2 1   3   1   2 1   2
diff:                                        (1,1)         (1,2)
change orient.:                               yes           no

                                              3 1           3 3
                                               *             *
                                              /|\           /|\
                                             / | \         / | \
                                            /  |  \       /  |  \
                                           /   |   \     /   |   \
                                          *----*----*   *----*----*
                                          1   2 2   3   1   2 2   1
diff:                                        (1,-1)        (1,1)
change orient.:                               no            yes

                                              3 1           3 2
                                               *             *
                                              /|\           /|\
                                             / | \         / | \
                                            /  |  \       /  |  \
                                           /   |   \     /   |   \
                                          *----*----*   *----*----*
                                          1   2 3   2   1   2 3   1
diff:                                        (1,-2)        (1,-1)
change orient.:                               yes           no


Change orientation if diff is:
(1,1) || (1,-2) || (2,2) || (2,-1) || (2,-1)

        */

        static void computeTriangleNeighbors (
            BoundaryGeometry* bg,
            std::vector<std::set<unsigned int>>& neighbors
            ) {
            std::vector<std::set<unsigned int>> adjacencies_to_pt (bg->Points_m.size());

            // for each triangles find adjacent triangles for each vertex
            for (unsigned int triangle_id = 0; triangle_id < bg->Triangles_m.size(); triangle_id++) {
                for (unsigned int j = 1; j <= 3; j++) {
                    auto pt_id = bg->PointID (triangle_id, j);
                    PAssert (pt_id < bg->Points_m.size ());
                    adjacencies_to_pt [pt_id].insert (triangle_id);
                }
            }

            for (unsigned int triangle_id = 0; triangle_id < bg->Triangles_m.size(); triangle_id++) {
                std::set<unsigned int>  to_A = adjacencies_to_pt [bg->PointID (triangle_id, 1)];
                std::set<unsigned int>  to_B = adjacencies_to_pt [bg->PointID (triangle_id, 2)];
                std::set<unsigned int>  to_C = adjacencies_to_pt [bg->PointID (triangle_id, 3)];

                std::set<unsigned int> intersect;
                std::set_intersection (
                    to_A.begin(), to_A.end(),
                    to_B.begin(), to_B.end(),
                    std::inserter(intersect,intersect.begin()));
                std::set_intersection(
                    to_B.begin(), to_B.end(),
                    to_C.begin(), to_C.end(),
                    std::inserter(intersect,intersect.begin()));
                std::set_intersection(
                    to_C.begin(), to_C.end(),
                    to_A.begin(), to_A.end(),
                    std::inserter(intersect, intersect.begin()));
                intersect.erase (triangle_id);

                neighbors [triangle_id] = intersect;
            }
            *gmsg << level2 << "* " << __func__ << ": Computing neighbors done" << endl;
        }

        /*
          Helper function for hasInwardPointingNormal()

          Determine if a point x is outside or inside the geometry or just on
          the boundary. Return true if point is inside geometry or on the
          boundary, false otherwise

          The basic idea here is:
          If a line segment from the point to test to a random point outside
          the geometry has has an even number of intersections with the
          boundary, the point is outside the geometry.

          Note:
          If the point is on the boundary, the number of intersections is 1.
          Points on the boundary are handled as inside.

          A random selection of the reference point outside the boundary avoids
          some specific issues, like line parallel to boundary.
         */
        static inline bool isInside (BoundaryGeometry* bg, const Vector_t x) {
            IpplTimings::startTimer (bg->TisInside_m);

            Vector_t y = Vector_t (
                bg->maxExtent_m[0] * (1.1 + gsl_rng_uniform(bg->randGen_m)),
                bg->maxExtent_m[1] * (1.1 + gsl_rng_uniform(bg->randGen_m)),
                bg->maxExtent_m[2] * (1.1 + gsl_rng_uniform(bg->randGen_m)));

            std::vector<Vector_t> intersection_points;
            //int num_intersections = 0;

            for (unsigned int triangle_id = 0; triangle_id < bg->Triangles_m.size(); triangle_id++) {
                Vector_t result;
                if (bg->intersectLineTriangle (SEGMENT, x, y, triangle_id, result)) {
                    intersection_points.push_back (result);
                    //num_intersections++;
                }
            }
            IpplTimings::stopTimer (bg->TisInside_m);
            return ((intersection_points.size () % 2) == 1);
        }

        // helper for function  makeTriangleNormalInwardPointing()
        static bool hasInwardPointingNormal (
            BoundaryGeometry* const bg,
            const int triangle_id
            ) {
            const Vector_t& A = bg->getPoint (triangle_id, 1);
            const Vector_t& B = bg->getPoint (triangle_id, 2);
            const Vector_t& C = bg->getPoint (triangle_id, 3);
            const Vector_t  triNormal = normalVector (A, B, C);

            // choose a point P close to the center of the triangle
            //const Vector_t P = (A+B+C)/3 + triNormal * 0.1;
            double minvoxelmesh = bg->voxelMesh_m.sizeOfVoxel[0];
            if (minvoxelmesh > bg->voxelMesh_m.sizeOfVoxel[1])
                minvoxelmesh = bg->voxelMesh_m.sizeOfVoxel[1];
            if (minvoxelmesh > bg->voxelMesh_m.sizeOfVoxel[2])
                minvoxelmesh = bg->voxelMesh_m.sizeOfVoxel[2];
            const Vector_t P = (A+B+C)/3 + triNormal * minvoxelmesh;
            /*
              The triangle normal points inward, if P is
              - outside the geometry and the dot product is negativ
              - or inside the geometry and the dot product is positiv

              Remember:
                The dot product is positiv only if both vectors are
                pointing in the same direction.
            */
            const bool is_inside = isInside (bg, P);
            const double dotPA_N = dot (P - A, triNormal);
            return (is_inside && dotPA_N >= 0) || (!is_inside && dotPA_N < 0);
        }

        // helper for function  makeTriangleNormalInwardPointing()
        static void orientTriangle (BoundaryGeometry* bg, int ref_id, int triangle_id) {
            // find pts of common edge
            int ic[2];
            int id[2];
            int n = 0;
            for (int i = 1; i <= 3; i++) {
                for (int j = 1; j <= 3; j++) {
                    if (bg->PointID (triangle_id, j) == bg->PointID (ref_id, i)) {
                        id[n] = j;
                        ic[n] = i;
                        n++;
                        if (n == 2) goto edge_found;
                    }
                }
            }
            PAssert (n == 2);
        edge_found:
            int diff = id[1] - id[0];
            if ((((ic[1] - ic[0]) == 1) && ((diff == 1) || (diff == -2))) ||
                (((ic[1] - ic[0]) == 2) && ((diff == -1) || (diff == 2)))) {
                std::swap (bg->PointID (triangle_id, id[0]), bg->PointID (triangle_id, id[1]));
            }
        }

        static void makeTriangleNormalInwardPointing (BoundaryGeometry* bg) {
            std::vector<std::set<unsigned int>> neighbors (bg->Triangles_m.size());

            computeTriangleNeighbors (bg, neighbors);

            // loop over all sub-meshes
            int triangle_id = 0;
            int parts = 0;
            std::vector<unsigned int> triangles (bg->Triangles_m.size());
            std::vector<unsigned int>::size_type queue_cursor = 0;
            std::vector<unsigned int>::size_type queue_end = 0;
            std::vector <bool> isOriented (bg->Triangles_m.size(), false);
            do {
                parts++;
                /*
                  Find next untested triangle, trivial for the first sub-mesh.
                  There is a least one not yet tested triangle!
                */
                while (isOriented[triangle_id])
                    triangle_id++;

                // ensure that normal of this triangle is inward pointing
                if (!hasInwardPointingNormal (bg, triangle_id)) {
                    std::swap (bg->PointID (triangle_id, 2), bg->PointID (triangle_id, 3));
                }
                isOriented[triangle_id] = true;

                // loop over all triangles in sub-mesh
                triangles[queue_end++] = triangle_id;
                do {
                    for (auto neighbor_id: neighbors[triangle_id]) {
                        if (isOriented[neighbor_id]) continue;
                        orientTriangle (bg, triangle_id, neighbor_id);
                        isOriented[neighbor_id] = true;
                        triangles[queue_end++] = neighbor_id;
                    }
                    queue_cursor++;
                } while (queue_cursor < queue_end && (triangle_id = triangles[queue_cursor],true));
            } while (queue_end < bg->Triangles_m.size());

            if (parts == 1) {
                *gmsg << level2 << "* " << __func__ << ": mesh is contiguous" << endl;
            } else {
                *gmsg << level2 << "* " << __func__ << ": mesh is discontiguous (" << parts << ") parts" << endl;
            }
            *gmsg << level2 <<"* Triangle Normal built done" << endl;
        }

    };

    debugFlags_m = 0;
    *gmsg << level2 << "* Initializing Boundary Geometry..." << endl;
    IpplTimings::startTimer (Tinitialize_m);

    if (!boost::filesystem::exists(h5FileName_m)) {
        throw OpalException("BoundaryGeometry::initialize",
                            "Failed to open file '" + h5FileName_m +
                            "', please check if it exists");
    }

    double xscale = Attributes::getReal(itsAttr[XSCALE]);
    double yscale = Attributes::getReal(itsAttr[YSCALE]);
    double zscale = Attributes::getReal(itsAttr[ZSCALE]);
    double xyzscale = Attributes::getReal(itsAttr[XYZSCALE]);
    double zshift = (double)(Attributes::getReal (itsAttr[ZSHIFT]));

    h5_int64_t rc;
#if defined (NDEBUG)
    (void)rc;
#endif
    rc = H5SetErrorHandler (H5AbortErrorhandler);
    PAssert (rc != H5_ERR);
    H5SetVerbosityLevel (1);

    h5_prop_t props = H5CreateFileProp ();
    MPI_Comm comm = Ippl::getComm();
    H5SetPropFileMPIOCollective (props, &comm);
    h5_file_t f = H5OpenFile (h5FileName_m.c_str(), H5_O_RDONLY, props);
    H5CloseProp (props);

    h5t_mesh_t* m = nullptr;
    H5FedOpenTriangleMesh (f, "0", &m);
    H5FedSetLevel (m, 0);

    auto numTriangles = H5FedGetNumElementsTotal (m);
    Triangles_m.resize (numTriangles);

    // iterate over all co-dim 0 entities, i.e. elements
    h5_loc_id_t local_id;
    int i = 0;
    h5t_iterator_t* iter = H5FedBeginTraverseEntities (m, 0);
    while ((local_id = H5FedTraverseEntities (iter)) >= 0) {
        h5_loc_id_t local_vids[4];
        H5FedGetVertexIndicesOfEntity (m, local_id, local_vids);
        PointID (i, 0) = 0;
        PointID (i, 1) = local_vids[0];
        PointID (i, 2) = local_vids[1];
        PointID (i, 3) = local_vids[2];
        i++;
    }
    H5FedEndTraverseEntities (iter);

    // loop over all vertices
    int num_points = H5FedGetNumVerticesTotal (m);
    Points_m.reserve (num_points);
    for (i = 0; i < num_points; i++) {
        h5_float64_t P[3];
        H5FedGetVertexCoordsByIndex (m, i, P);
        Points_m.push_back (Vector_t (
            P[0] * xyzscale * xscale,
            P[1] * xyzscale * yscale,
            P[2] * xyzscale * zscale + zshift));
    }
    H5FedCloseMesh (m);
    H5CloseFile (f);
    *gmsg << level2 << "* Reading mesh done" << endl;

    Local::computeGeometryInterval (this);
    computeMeshVoxelization ();
    haveInsidePoint_m = false;
    std::vector<double> pt = Attributes::getRealArray (itsAttr[INSIDEPOINT]);
    if (!pt.empty()) {
        if (pt.size () != 3) {
            throw OpalException (
                "BoundaryGeometry::initialize()",
                "Dimension of INSIDEPOINT must be 3");
        }
        /* test whether this point is inside */
        insidePoint_m = {pt[0], pt[1], pt[2]};
        bool is_inside = isInside (insidePoint_m);
        if (is_inside == false) {
            throw OpalException (
                "BoundaryGeometry::initialize()",
                "INSIDEPOINT is not inside the geometry");
        }
        haveInsidePoint_m = true;
    } else {
        haveInsidePoint_m = findInsidePoint();
    }
    if (haveInsidePoint_m == true) {
        *gmsg << level2 << "* using as point inside the geometry: ("
              << insidePoint_m[0] << ", "
              << insidePoint_m[1] << ", "
              << insidePoint_m[2] << ")" << endl;
    } else {
        *gmsg << level2 << "* no point inside the geometry found!" << endl;
    }

    Local::makeTriangleNormalInwardPointing (this);

    TriNormals_m.resize (Triangles_m.size());
    TriAreas_m.resize (Triangles_m.size());

    for (size_t i = 0; i < Triangles_m.size(); i++) {
        const Vector_t& A = getPoint (i, 1);
        const Vector_t& B = getPoint (i, 2);
        const Vector_t& C = getPoint (i, 3);

        TriAreas_m[i] = computeArea (A, B, C);
        TriNormals_m[i] = normalVector (A, B, C);

    }
    *gmsg << level2 << "* Triangle barycent built done" << endl;

    *gmsg << *this << endl;
    Ippl::Comm->barrier();
    IpplTimings::stopTimer (Tinitialize_m);
}

/*
  Line segment triangle intersection test. This method should be used only
  for "tiny" line segments or, to be more exact, if the number of
  voxels covering the bounding box of the line segment is small (<<100).

  Actually the method can be used for any line segment, but may not perform
  well. Performace depends on the size of the bounding box of the line
  segment.

  The method returns the number of intersections of the line segment defined
  by the points P and Q with the boundary. If there are multiple intersections,
  the nearest intersection point with respect to P wil be returned.
 */
int
BoundaryGeometry::intersectTinyLineSegmentBoundary (
    const Vector_t& P,                  // [i] starting point of ray
    const Vector_t& Q,                  // [i] end point of ray
    Vector_t& intersect_pt,             // [o] intersection with boundary
    int& triangle_id                    // [o] intersected triangle
    ) {
#ifdef ENABLE_DEBUG
    if (debugFlags_m & debug_intersectTinyLineSegmentBoundary) {
        *gmsg << "* " << __func__ << ": "
              << "  P = " << P
              << ", Q = " << Q
              << endl;
    }
#endif
    const Vector_t v_ = Q - P;
    const Ray r = Ray (P, v_);
    const Vector_t bbox_min = {
        std::min(P[0], Q[0]),
        std::min(P[1], Q[1]),
        std::min(P[2], Q[2]) };
    const Vector_t bbox_max = {
        std::max(P[0], Q[0]),
        std::max(P[1], Q[1]),
        std::max(P[2], Q[2]) };
    int i_min, i_max;
    int j_min, j_max;
    int k_min, k_max;
    mapPoint2VoxelIndices (bbox_min, i_min, j_min, k_min);
    mapPoint2VoxelIndices (bbox_max, i_max, j_max, k_max);

    Vector_t tmp_intersect_pt = Q;
    double tmin = 1.1;

    /*
      Triangles can - and in many cases do - intersect with more than one
      voxel.  If we loop over all voxels intersecting with the line segment
      spaned by the points P and Q, we might perform the same line-triangle
      intersection test more than once.  We must this into account when
      counting the intersections with the boundary.

      To avoid multiple counting we can either
      - build a set of all relevant triangles and loop over this set
      - or we loop over all voxels and remember the intersecting triangles.

      The first solution is implemented here.
     */
    std::unordered_set<int> triangle_ids;
    for (int i = i_min; i <= i_max; i++) {
        for (int j = j_min; j <= j_max; j++) {
            for (int k = k_min; k <= k_max; k++) {
                const Vector_t bmin = mapIndices2Voxel(i, j, k);
                const Voxel v(bmin, bmin + voxelMesh_m.sizeOfVoxel);
#ifdef ENABLE_DEBUG
                if (debugFlags_m & debug_intersectTinyLineSegmentBoundary) {
                    *gmsg << "* " << __func__ << ": "
                          << "  Test voxel: (" << i << ", " << j << ", " << k << "), "
                          << v.pts[0] << v.pts[1]
                          << endl;
                }
#endif
                /*
                  do line segment and voxel intersect? continue if not
                */
                if (!v.intersect (r)) {
                    continue;
                }

                /*
                  get triangles intersecting with this voxel and add them to
                  the to be tested triangles.
                 */
                const int voxel_id = mapVoxelIndices2ID (i, j, k);
                const auto triangles_intersecting_voxel =
                    voxelMesh_m.ids.find (voxel_id);
                if (triangles_intersecting_voxel != voxelMesh_m.ids.end()) {
                    triangle_ids.insert (
                        triangles_intersecting_voxel->second.begin(),
                        triangles_intersecting_voxel->second.end());
                }
            }
        }
    }
    /*
      test all triangles intersecting with one of the above voxels
      if there is more than one intersection, return closest
    */
    int num_intersections = 0;
    int tmp_intersect_result = 0;

    for (auto it = triangle_ids.begin ();
         it != triangle_ids.end ();
         it++) {

        tmp_intersect_result = intersectLineTriangle (
            LINE,
            P, Q,
            *it,
            tmp_intersect_pt);
#ifdef ENABLE_DEBUG
        if (debugFlags_m & debug_intersectTinyLineSegmentBoundary) {
            *gmsg << "* " << __func__ << ": "
                  << "  Test triangle: " << *it
                  << "  intersect: " << tmp_intersect_result
                  << getPoint(*it,1)
                  << getPoint(*it,2)
                  << getPoint(*it,3)
                  << endl;
        }
#endif
        switch (tmp_intersect_result) {
        case 0:                     // no intersection
        case 2:                     // both points are outside
        case 4:                     // both points are inside
            break;
        case 1:                     // line and triangle are in same plane
        case 3:                     // unique intersection in segment
            double t;
            if (cmp::eq_zero(Q[0] - P[0]) == false) {
                t = (tmp_intersect_pt[0] - P[0]) / (Q[0] - P[0]);
            } else if (cmp::eq_zero(Q[1] - P[1]) == false) {
                t = (tmp_intersect_pt[1] - P[1]) / (Q[1] - P[1]);
            } else {
                t = (tmp_intersect_pt[2] - P[2]) / (Q[2] - P[2]);
            }
            num_intersections++;
            if (t < tmin) {
#ifdef ENABLE_DEBUG
                if (debugFlags_m & debug_intersectTinyLineSegmentBoundary) {
                    *gmsg << "* " << __func__ << ": "
                          << "  set triangle"
                          << endl;
                }
#endif
                tmin = t;
                intersect_pt = tmp_intersect_pt;
                triangle_id = (*it);
            }
            break;
        case -1:                    // triangle is degenerated
            PAssert (tmp_intersect_result != -1);
            exit (42);              // terminate even if NDEBUG is set
        }
    }                   // end for all triangles
    return num_intersections;
}

/*
  General purpose line segment boundary intersection test.

  The method returns with a value > 0 if an intersection was found.
 */
int
BoundaryGeometry::intersectLineSegmentBoundary (
    const Vector_t& P0,                 // [in] starting point of ray
    const Vector_t& P1,                 // [in] end point of ray
    Vector_t& intersect_pt,             // [out] intersection with boundary
    int& triangle_id                    // [out] triangle the line segment intersects with
    ) {
#ifdef ENABLE_DEBUG
    int saved_flags = debugFlags_m;
    if (debugFlags_m & debug_intersectLineSegmentBoundary) {
        *gmsg << "* " << __func__ << ": "
              << "  P0 = " << P0
              << "  P1 = " << P1
              << endl;
        debugFlags_m |= debug_intersectTinyLineSegmentBoundary;
    }
#endif
    triangle_id = -1;

    const Vector_t v = P1 - P0;
    int intersect_result = 0;
    int n = 0;
    int i_min, j_min, k_min;
    int i_max, j_max, k_max;
    do {
        n++;
        Vector_t Q = P0 + v / n;
        Vector_t bbox_min = {
            std::min(P0[0], Q[0]),
            std::min(P0[1], Q[1]),
            std::min(P0[2], Q[2]) };
        Vector_t bbox_max = {
            std::max(P0[0], Q[0]),
            std::max(P0[1], Q[1]),
            std::max(P0[2], Q[2]) };
        mapPoint2VoxelIndices (bbox_min, i_min, j_min, k_min);
        mapPoint2VoxelIndices (bbox_max, i_max, j_max, k_max);
    } while (( (i_max-i_min+1) * (j_max-j_min+1) * (k_max-k_min+1)) > 27);
    Vector_t P = P0;
    Vector_t Q;
    const Vector_t v_ = v / n;

    for (int l = 1; l <= n; l++, P = Q) {
        Q = P0 + l*v_;
        intersect_result = intersectTinyLineSegmentBoundary (
            P, Q, intersect_pt, triangle_id);
        if (triangle_id != -1) {
            break;
        }
    }
#ifdef ENABLE_DEBUG
    if (debugFlags_m & debug_intersectLineSegmentBoundary) {
        *gmsg << "* " << __func__ << ": "
              << "  result=" << intersect_result
              << "  intersection pt: " << intersect_pt
              << endl;
        debugFlags_m = saved_flags;
    }
#endif
    return intersect_result;
}

/**
   Determine whether a particle with position @param r, momenta @param v,
   and time step @param dt will hit the boundary.

   return value:
        -1  no collison with boundary
        0   particle will collide with boundary in next time step
 */
int
BoundaryGeometry::partInside (
    const Vector_t& r,                  // [in] particle position
    const Vector_t& v,                  // [in] momentum
    const double dt,                    // [in]
    Vector_t& intersect_pt,             // [out] intersection with boundary
    int& triangle_id                    // [out] intersected triangle
    ) {
#ifdef ENABLE_DEBUG
    int saved_flags = debugFlags_m;
    if (debugFlags_m & debug_PartInside) {
        *gmsg << "* " << __func__ << ": "
              << "  r=" << r
              << "  v=" << v
              << "  dt=" << dt
              << endl;
        debugFlags_m |= debug_intersectTinyLineSegmentBoundary;
    }
#endif
    int ret = -1;                       // result defaults to no collision

    // nothing to do if momenta == 0
    if (v == (Vector_t)0)
        return ret;

    IpplTimings::startTimer (TPartInside_m);

    // P0, P1: particle position in time steps n and n+1
    const Vector_t P0 = r;
    const Vector_t P1 = r + (Physics::c * v * dt / std::sqrt (1.0 + dot(v,v)));

    Vector_t tmp_intersect_pt = 0.0;
    int tmp_triangle_id = -1;
    intersectTinyLineSegmentBoundary (P0, P1, tmp_intersect_pt, tmp_triangle_id);
    if (tmp_triangle_id >= 0) {
        intersect_pt = tmp_intersect_pt;
        triangle_id = tmp_triangle_id;
        ret = 0;
    }
#ifdef ENABLE_DEBUG
    if (debugFlags_m & debug_PartInside) {
        *gmsg << "* " << __func__ << ":"
              << "  result=" << ret;
        if (ret == 0) {
            *gmsg << "  intersetion=" << intersect_pt;
        }
        *gmsg << endl;
        debugFlags_m = saved_flags;
    }
#endif
    IpplTimings::stopTimer (TPartInside_m);
    return ret;
}

void
BoundaryGeometry::writeGeomToVtk (std::string fn) {
    std::ofstream of;
    of.open (fn.c_str ());
    PAssert (of.is_open ());
    of.precision (6);
    of << "# vtk DataFile Version 2.0" << std::endl;
    of << "generated using DataSink::writeGeoToVtk" << std::endl;
    of << "ASCII" << std::endl << std::endl;
    of << "DATASET UNSTRUCTURED_GRID" << std::endl;
    of << "POINTS " << Points_m.size () << " float" << std::endl;
    for (unsigned int i = 0; i < Points_m.size (); i++)
        of << Points_m[i](0) << " "
           << Points_m[i](1) << " "
           << Points_m[i](2) << std::endl;
    of << std::endl;

    of << "CELLS "
       << Triangles_m.size() << " "
       << 4 * Triangles_m.size() << std::endl;
    for (size_t i = 0; i < Triangles_m.size(); i++)
        of << "3 "
           << PointID (i, 1) << " "
           << PointID (i, 2) << " "
           << PointID (i, 3) << std::endl;
    of << "CELL_TYPES " << Triangles_m.size() << std::endl;
    for (size_t i = 0; i < Triangles_m.size(); i++)
        of << "5" << std::endl;
    of << "CELL_DATA " << Triangles_m.size() << std::endl;
    of << "SCALARS " << "cell_attribute_data" << " float " << "1" << std::endl;
    of << "LOOKUP_TABLE " << "default" << std::endl;
    for (size_t i = 0; i < Triangles_m.size(); i++)
        of << (float)(i) << std::endl;
    of << std::endl;
}

Inform&
BoundaryGeometry::printInfo (Inform& os) const {
    os << endl;
    os << "* ************* B O U N D A R Y  G E O M E T R Y *********************************** " << endl;
    os << "* GEOMETRY                  "  << getOpalName () << '\n'
       << "* FGEOM                     '" << Attributes::getString (itsAttr[FGEOM])  << "'\n"
       << "* TOPO                      " << Attributes::getString (itsAttr[TOPO])   << '\n'
       << "* XSCALE                    " << Attributes::getReal (itsAttr[XSCALE])   << '\n'
       << "* YSCALE                    " << Attributes::getReal (itsAttr[YSCALE])   << '\n'
       << "* ZSCALE                    " << Attributes::getReal (itsAttr[ZSCALE])   << '\n'
       << "* XYZSCALE                  " << Attributes::getReal (itsAttr[XYZSCALE]) << '\n'
       << "* LENGTH                    " << Attributes::getReal (itsAttr[LENGTH])   << '\n'
       << "* S                         " << Attributes::getReal (itsAttr[S]) << '\n'
       << "* A                         " << Attributes::getReal (itsAttr[A]) << '\n'
       << "* B                         " << Attributes::getReal (itsAttr[B]) << '\n';
    if (getTopology () == Topology::BOXCORNER) {
        os << "* C                         " << Attributes::getReal (itsAttr[C]) << '\n'
           << "* L1                        " << Attributes::getReal (itsAttr[L1]) << '\n'
           << "* L2                        " << Attributes::getReal (itsAttr[L2]) << '\n';
    }
    os << "* Total triangle num        " << Triangles_m.size() << '\n'
       << "* Total points num          " << Points_m.size () << '\n'
       << "* Geometry bounds(m) Max =  " << maxExtent_m << '\n'
       << "*                    Min =  " << minExtent_m << '\n'
       << "* Geometry length(m)        " << maxExtent_m - minExtent_m << '\n'
       << "* Resolution of voxel mesh  " << voxelMesh_m.nr_m << '\n'
       << "* Size of voxel             " << voxelMesh_m.sizeOfVoxel << '\n'
       << "* Number of voxels in mesh  " << voxelMesh_m.ids.size () << endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}
