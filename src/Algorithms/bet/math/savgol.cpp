/*
   Project: Beam Envelope Tracker (BET)
   Author:  Rene Bakker et al.
   Created: 23-11-2006

   Savitzky-Golay Smoothing Filters

   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43108-5)
*/

#include <algorithm>
#include <cmath>
#include <cstring>

#include "Utilities/OpalException.h"
#include "Algorithms/bet/math/savgol.h"

// see: https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
        // unless the result is subnormal
        || std::abs(x-y) < std::numeric_limits<T>::min();
}

/* internal functions
   ====================================================================== */

/* lubksb()
   Solves the set of n linear equations A·X = B.
   Here a[1..n][1..n] is input, not as the matrix A but rather as its
   LU decomposition, determined by the routine ludcmp. indx[1..n] is input
   as the permutation vector returned by ludcmp. b[1..n] is input as the
   right-hand side vector B, and returns with the solution vector X. a, n,
   and indx are not modified by this routine and can be left in place for
   successive calls with different right-hand sides b. This routine takes
   into account the possibility that b will begin with many zero elements,
   so it is efficient for use in matrix inversion.
*/
static void lubksb(double **a, int n, int indx[], double b[]) {

    for(auto i = 1; i <= n; i++) {
        auto ip = indx[i];
        auto sum = b[ip];
        auto ii = 0;
        b[ip] = b[i];
        if(ii)
            for(auto j = ii; j <= i - 1; j++)
                sum -= a[i][j] * b[j];
        else if(sum)
            ii = i;
        b[i] = sum;
    }
    for(auto i = n; i >= 1; i--) {
        auto sum = b[i];
        for(auto j = i + 1; j <= n; j++)
            sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}

#define TINY 1.0e-20;

/* ludcmp()
   Given a matrix a[1..n][1..n], this routine replaces it by the LU
   decomposition of a rowwise permutation of itself. a and n are input.
   a is output, arranged as in equation (2.3.14) above;
   indx[1..n] is an output vector that records the row permutation effected
   by the partial pivoting; d is output as ±1 depending on whether the number
   of row interchanges was even or odd, respectively. This routine is used
   in combination with lubksb to solve linear equations or invert a matrix.
*/
static void ludcmp(double **a, int n, int indx[], double &d) {

    double vv[n];
    d = 1.0;
    for(auto i = 1; i <= n; i++) {
        double big = 0.0;
        for(auto j = 1; j <= n; j++) {
            auto temp = std::abs(a[i][j]);
            if(temp > big)
                big = temp;
        }
        if (!almost_equal (big, 0.0, 4)) {
            throw OpalException (
                "BetMath_savgol()",
                "Singular matrix.");
        }
        vv[i] = 1.0 / big;
    }
    for(auto j = 1; j <= n; j++) {
        for(auto i = 1; i < j; i++) {
            auto sum = a[i][j];
            for(auto k = 1; k < i; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        double big = 0.0;
        int imax = -1;
        for(auto i = j; i <= n; i++) {
            auto sum = a[i][j];
            for(auto k = 1; k < j; k++)
                sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            auto temp = vv[i] * std::abs(sum);
            if(temp >= big) {
                big = temp;
                imax = i;
            }
        }
        if(j != imax) {
            for(auto k = 1; k <= n; k++) {
                auto temp = a[imax][k];
                a[imax][k] = a[j][k];
                a[j][k] = temp;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if(a[j][j] <= std::numeric_limits<double>::epsilon())
            a[j][j] = TINY;
        if(j != n) {
            double temp = 1.0 / (a[j][j]);
            for(auto i = j + 1; i <= n; i++)
                a[i][j] *= temp;
        }
    }
}

/* savgol()
   Returns in c[1..np], in wrap-around order (N.B.!)
   consistent with the argument respons in
   routine convlv, a set of Savitzky-Golay filter coefficients.
   nl is the number of leftward (past) data points used, while
   nr is the number of rightward (future) data points, making
      the total number of data points used nl+nr+1.
   ld is the order of the derivative desired
      (e.g., ld = 0 for smoothed function).
   m  is the order of the smoothing polynomial,
      also equal to the highest conserved moment;
      usual values are m = 2 or m = 4.
*/
static void savgol(double c[], int np, int nl, int nr, int ld, int m) {

    if (np >= nl + nr + 1 && nl >= 0 && nr >= 0 && ld <= m && nl + nr >= m) {
        throw OpalException(
            "BetMath_savgol()",
            "Bad arguments.");
    }
    int indx[m + 1];

    // allocate matrix on the stack. Since the matrix is small, this is ok.
    auto nrow = m + 2;
    auto ncol = nrow;
    double* a[nrow+1];
    double matrix[nrow * ncol + 1];
    a[0] = matrix;
    for (auto i = 1; i <= m+1; i++) {
         a[i] = a[i-1] + ncol;
    }

    double b[m + 1];
    for(int ipj = 0; ipj <= (m << 1); ipj++) {
        auto sum = (ipj ? 0.0 : 1.0);
        for(auto k = 1; k <= nr; k++)
            sum += std::pow((double)k, ipj);
        for(auto k = 1; k <= nl; k++)
            sum += std::pow((double) -k, ipj);
        auto mm = std::min(ipj, 2 * m - ipj);
        for(auto imj = -mm; imj <= mm; imj += 2)
            a[1+(ipj+imj)/2][1+(ipj-imj)/2] = sum;
    }
    double d = 1.0;
    ludcmp(a, m + 1, indx, d);
    for(auto j = 1; j <= m + 1; j++)
        b[j] = 0.0;
    b[ld+1] = 1.0;
    lubksb(a, m + 1, indx, b);
    for(auto kk = 1; kk <= np; kk++)
        c[kk] = 0.0;
    for(auto k = -nl; k <= nr; k++) {
        auto sum = b[1];
        double fac = 1.0;
        for(auto mm = 1; mm <= m; mm++)
            sum += b[mm+1] * (fac *= k);
        auto kk = ((np - k) % np) + 1;
        c[kk] = sum;
    }
}

/* four1()
   Replaces data[1..2*nn] by its discrete Fourier transform, if isign
   is input as 1; or replaces data[1..2*nn] by nn times its inverse
   discrete Fourier transform, if isign is input as -1. data is a
   complex array of length nn or, equivalently, a real array of length
   2*nn. nn MUST be an integer power of 2 (this is not checked for!).
*/
static void four1(double data[], int nn, int isign) {

    auto n = nn << 1;
    auto j = 1;
    for(auto i = 1; i < n; i += 2) {
        if(j > i) {
            std::swap(data[j], data[i]);
            std::swap(data[j+1], data[i+1]);
        }
        auto m = n >> 1;
        while(m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    auto mmax = 2;
    while(n > mmax) {
        auto istep = mmax << 1;
        auto theta = isign * (6.28318530717959 / mmax);
        auto wtemp = sin(0.5 * theta);
        auto wpr = -2.0 * wtemp * wtemp;
        auto wpi = sin(theta);
        auto wr = 1.0;
        auto wi = 0.0;
        for(auto m = 1; m < mmax; m += 2) {
            for(auto i = m; i <= n; i += istep) {
                j = i + mmax;
                auto tempr = wr * data[j] - wi * data[j+1];
                auto tempi = wr * data[j+1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j+1] = data[i+1] - tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}

/* realft()
   Calculates the Fourier transform of a set of n real-valued data points.
   Replaces this data (which is stored in array data[1..n]) by the positive
   frequency half of its complex Fourier transform. The real-valued first
   and last components of the complex transform are returned as elements
   data[1] and data[2], respectively. n must be a power of 2. This routine
   also calculates the inverse transform of a complex data array if it is
   the transform of real data. (Result in this case must be multiplied
   by 2/n.)static
*/
static void realft(double data[], int n, int isign) {
    const double c1 = 0.5;
    double c2;

    auto theta = 3.141592653589793 / (double)(n >> 1);
    if(isign == 1) {
        c2 = -0.5;
        four1(data, n >> 1, 1);
    } else {
        c2 = 0.5;
        theta = -theta;
    }
    auto wtemp = sin(0.5 * theta);
    auto wpr = -2.0 * wtemp * wtemp;
    auto wpi = sin(theta);
    auto wr = 1.0 + wpr;
    auto wi = wpi;
    auto np3 = n + 3;
    double h1r;
    for(auto i = 2; i <= (n >> 2); i++) {
        auto i1 = 2 * i - 1;
        auto i2 = 2 * i;
        auto i3 = np3 - i2;
        auto i4 = 1 + i3;
        h1r = c1 * (data[i1] + data[i3]);
        auto h1i = c1 * (data[i2] - data[i4]);
        auto h2r = -c2 * (data[i2] + data[i4]);
        auto h2i = c2 * (data[i1] - data[i3]);
        data[i1] = h1r + wr * h2r - wi * h2i;
        data[i2] = h1i + wr * h2i + wi * h2r;
        data[i3] = h1r - wr * h2r + wi * h2i;
        data[i4] = -h1i + wr * h2i + wi * h2r;
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    if(isign == 1) {
        data[1] = (h1r = data[1]) + data[2];
        data[2] = h1r - data[2];
    } else {
        data[1] = c1 * ((h1r = data[1]) + data[2]);
        data[2] = c1 * (h1r - data[2]);
        four1(data, n >> 1, -1);
    }
}

/* twofft()
   Given two real input arrays data1[1..n] and data2[1..n],
   this routine calls four1 and returns two complex output arrays,
   fft1[1..2n] and fft2[1..2n], each of complex length n (i.e.,
   real length 2*n), which contain the discrete Fourier transforms
   of the respective data arrays. n MUST be an integer power of 2
*/
static void twofft(double data1[], double data2[], double fft1[], double fft2[],
                   int n) {

    auto nn2 = 2 + n + n;
    auto nn3 = 1 + nn2;
    for(auto j = 1, jj = 2; j <= n; j++, jj += 2) {
        fft1[jj-1] = data1[j];
        fft1[jj] = data2[j];
    }
    four1(fft1, n, 1);
    fft2[1] = fft1[2];
    fft1[2] = fft2[2] = 0.0;
    for(auto j = 3; j <= n + 1; j += 2) {
        auto rep = 0.5 * (fft1[j] + fft1[nn2-j]);
        auto rem = 0.5 * (fft1[j] - fft1[nn2-j]);
        auto aip = 0.5 * (fft1[j+1] + fft1[nn3-j]);
        auto aim = 0.5 * (fft1[j+1] - fft1[nn3-j]);
        fft1[j] = rep;
        fft1[j+1] = aim;
        fft1[nn2-j] = rep;
        fft1[nn3-j] = -aim;
        fft2[j] = aip;
        fft2[j+1] = -rem;
        fft2[nn2-j] = aip;
        fft2[nn3-j] = rem;
    }
}

/* convlv()
   Convolves or deconvolves a real data set data[1..n] (including any
   user-supplied zero padding) with a response function respns[1..n].
   The response function must be stored in wrap-around order in the
   first m elements of respns, where m is an odd integer <=n. Wrap-
   around order means that the first half of the array respns contains
   the impulse response function at positive times, while the second
   half of the array contains the impulse response function at negative
   times, counting down from the highest element respns[m]. On input
   isign is +1 for convolution, -1 for deconvolution. The answer is
   returned in the first n components of ans. However, ans must be
   supplied in the calling program with dimensions [1..2*n], for
   consistency with twofft. n MUST be an integer power of two.
*/
static void convlv(double data[], int n, double respns[], int m,
            int isign, double ans[]) {

    double fft[n << 1];

    for(auto i = 1; i <= (m - 1) / 2; i++)
        respns[n+1-i] = respns[m+1-i];
    for(auto i = (m + 3) / 2; i <= n - (m - 1) / 2; i++)
        respns[i] = 0.0;
    twofft(data, respns, fft, ans, n);
    auto no2 = n >> 1;
    for(auto i = 2; i <= n + 2; i += 2) {
        auto dum = ans[i-1];
        if(isign == 1) {
            ans[i-1] = (fft[i-1] * dum - fft[i] * ans[i]) / no2;
            ans[i] = (fft[i] * dum + fft[i-1] * ans[i]) / no2;
        } else {
            auto mag2 = std::pow(ans[i-1], 2) + std::pow(ans[i], 2);
            if(!almost_equal (mag2, 0.0, 4)) {
                throw OpalException(
                    "BetMath_savgol::convlv()",
                    "Deconvoloving at repsonse zero.");
            }
            ans[i-1] = (fft[i-1] * dum + fft[i] * ans[i]) / mag2 / no2;
            ans[i] = (fft[i] * dum - fft[i-1] * ans[i]) / mag2 / no2;
        }
    }
    ans[2] = ans[n+1];
    realft(ans, n, -1);
}

/* sgSmooth()
   Smoothes c[0..n-1] with a Savitzky-Golay filter.
   nl is the number of leftward (past) data points used, while
   nr is the number of rightward (future) data points, making
      the total number of data points used nl+nr+1.
   ld is the order of the derivative desired
      (e.g., ld = 0 for smoothed function).
   m  is the order of the smoothing polynomial,
      also equal to the highest conserved moment;
      usual values are m = 2 or m = 4.
*/
void sgSmooth(double c[], int n, int nl, int nr, int ld, int m) {

    // make dimension 2^m with m integer
    auto temp = n;
    auto log2_of_n = 0;
    while (temp >>= 1) ++log2_of_n;
    auto nn = int(std::pow(2, log2_of_n + 1));

    // memory allocation
    double cf[nn];
    double cIn[nn];
    double cOut[nn * 2];

    // fill data array
    cIn[0] = 0.0;
    memcpy(&cIn[1], c, sizeof(double)*n);
    for(auto i = n + 1; i <= nn; i++) {
        cIn[i] = 0.0;
    }

    // create filter coefficients
    auto np = nl + nr + 1;
    if((np % 2) == 0) ++np;
    savgol(cf, np, nl, nr, ld, m);

    // filter
    int isign = 1;
    convlv(cIn, nn, cf, np, isign, cOut);

    // move data back
    memcpy(c, &cOut[1], sizeof(double)*n);
}

// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
