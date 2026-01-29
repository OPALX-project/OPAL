#ifndef CMPFLOAT_H
#define CMPFLOAT_H

#include <cfloat>
#include <cmath>

/*
  In the following namespaces various approximately floating point
  comparisons are implemented. The used implementation is selected
  via

  namespaces cmp = IMPLEMENTATION;

*/

/*
  macros used in all implementations/namespaces
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

#endif
