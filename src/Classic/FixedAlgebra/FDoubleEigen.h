#ifndef CLASSIC_FDoubleEigen_HH
#define CLASSIC_FDoubleEigen_HH

// ------------------------------------------------------------------------
// $RCSfile: FDoubleEigen.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: FDoubleEigen
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:36 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include <complex>

// Class FDoubleEigen
// ------------------------------------------------------------------------
/// Eigenvalues and eigenvectors for a real general matrix.
//  Translated to FORTRAN by Burton S. Garbow, ANL.
//  Adapted March 1997 by F. Christoph Iselin, CERN, SL/AP.
//  Changed December 1998 to template by F. Christoph Iselin, CERN, SL/AP.

template <int N>
class FDoubleEigen {

public:

    /// Constructor.
    //  Find eigenvalues for matrix [b]M[/b].
    //  If [b]vec[/b] is true, the eigenvectors are also computed.
    FDoubleEigen(const FMatrix<double, N, N> &M, bool vec = false);

    FDoubleEigen();
    FDoubleEigen(const FDoubleEigen &);
    ~FDoubleEigen();

    /// Get eigenvalues.
    //  Return eigenvalues as a complex vector.
    FVector<std::complex<double>, N> eigenValues() const;

    /// Get eigenvectors.
    //  Return eigenvectors as a complex matrix.
    FMatrix<std::complex<double>, N, N> eigenVectors() const;

    /// Get eigenvectors.
    //  Return eigenvectors packed in a real matrix.
    //  Real eigenvectors appear a single column in the same position as
    //  the corresponding real eigenvalue.
    //  Complex eigenvalues occur in conjugate pairs, and the corresponding
    //  eigenvectors appear as two real columns containing the real and
    //  imaginary parts respectively.
    FMatrix<double, N, N> packedEigenVectors() const;

private:

    // Not implemented.
    void operator=(const FDoubleEigen &);

    // Used by eigenvalue and eigenvector routines
    static void balance(FMatrix<double, N, N> &, int &low, int &high,
                        double scale[N]);

    static void exchange(FMatrix<double, N, N> &, int j, int m,
                         int low, int high);

    static void elmhes(FMatrix<double, N, N> &, int low,
                       int high, int index[N]);

    void elmtran(FMatrix<double, N, N> &, int low, int high,
                 int index[N]);

    int hqr(FMatrix<double, N, N> &, int low, int high);

    int hqr2(FMatrix<double, N, N> &, int low, int high);

    void balbak(int low, int high, double scale[N]);

    // Representation of the eigenvalues and eigenvectors.
    FVector<std::complex<double>, N> lambda;
    FMatrix<double, N, N>       vectors;
};


// Implementation.
// ------------------------------------------------------------------------

#include "Utilities/EigenvalueError.h"
#include <algorithm>
#include <cmath>

inline
void cdiv(double ar, double ai, double br, double bi,
          double &cr, double &ci)
// Complex division: (cr,ci) = (ar,ai) / (br,bi).
// Adapted March 1997 by F. Christoph Iselin, CERN, SL/AP.
{
    double s = std::abs(br) + std::abs(bi);
    double ars = ar / s;
    double ais = ai / s;
    double brs = br / s;
    double bis = bi / s;
    
    s = brs * brs + bis * bis;
    
    cr = (ars * brs + ais * bis) / s;
    ci = (ais * brs - ars * bis) / s;
}

// Class FDoubleEigen; public methods.
// --------------------------------------------------------------------------

template <int N>
FDoubleEigen<N>::FDoubleEigen():
    lambda(), vectors()
{}


template <int N>
FDoubleEigen<N>::FDoubleEigen(const FDoubleEigen &M):
    lambda(M.lambda), vectors(M.vectors)
{}


template <int N>
FDoubleEigen<N>::FDoubleEigen(const FMatrix<double, N, N> &M, bool vec):
    lambda(), vectors()
    // Before execution:
    //   M        is the matrix whose eigenvalues are requested.
    // After execution:
    //   lambda   contains the eigenvalues.
    // If vect is true, then:
    //   vectors  contains the eigenvectors, packed in a real matrix.
{
    for(int i = 0; i < N; ++i) vectors(i, i) = 1.0;
    int low, upp;

    // Copy original matrix.
    FMatrix<double, N, N> h(M);

    // Balance the copy.
    double scale[N];
    balance(h, low, upp, scale);

    // Reduce the copy to upper Hessenberg form.
    int index[N];
    elmhes(h, low, upp, index);

    if(vec) {
        elmtran(h, low, upp, index);

        if(hqr2(h, low, upp) != 0) {
            throw EigenvalueError("FDoubleEigen::FDoubleEigen()",
                                  "Unable to find all eigenvalues.");
        } else {
            balbak(low, upp, scale);
        }
    } else {
        if(hqr(h, low, upp) != 0) {
            throw EigenvalueError("FDoubleEigen::FDoubleEigen()",
                                  "Unable to find all eigenvalues.");
        }
    }
}


template <int N>
FDoubleEigen<N>::~FDoubleEigen()
{}


template <int N>
FVector<std::complex<double>, N> FDoubleEigen<N>::eigenValues() const
// Return Eigenvalues as complex vector.
{
    return lambda;
}


template <int N>
FMatrix<std::complex<double>, N, N> FDoubleEigen<N>::eigenVectors() const
// Return eigenvectors as a complex matrix.
{
    FMatrix<std::complex<double>, N, N> R;

    for(int i = 0; i < N; i++) {
        if(std::imag(lambda[i]) == 0.0) {
            // One real eigenvector.
            for(int j = 0; j < N; j++) {
                R[j][i] = std::complex<double>(vectors[j][i]);
            }
        } else {
            // Two complex eigenvectors.
            for(int j = 0; j < N; j++) {
                R[j][i]   = std::complex<double>(vectors[j][i], + vectors[j][i+1]);
                R[j][i+1] = std::complex<double>(vectors[j][i], - vectors[j][i+1]);
            }
            i++;
        }
    }

    return R;
}


template <int N>
FMatrix<double, N, N> FDoubleEigen<N>::packedEigenVectors() const
// Return eigenvectors packed in a real matrix.
{
    return vectors;
}


// Class FDoubleEigen; private methods.
// --------------------------------------------------------------------------

template <int N>
void FDoubleEigen<N>::balance
(FMatrix<double, N, N> &copy, int &low, int &upp, double scale[N])

// This method is a translation of the Algol procedure balance,
// Num. Math. 13, 293-304(1969) by Parlett and Reinsch.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 315-326(1971).
//
// It balances a real matrix and isolates eigenvalues whenever possible.
//
// On input:
// copy       Contains a copy of the original matrix A.
//
// On output:
// copy       Contains the balanced matrix A.
//
// low, upp   Are two indices such that A[i][j] is equal to zero if
//            (1) i is greater than j and
//            (2) j = 0, ..., low - 1 or i = upp + 1, ..., N - 1.
//
// scale      Contains information determining the
//            permutations and scaling factors used.
//
// Suppose that the principal submatrix in rows "low" through "upp"
// has been balanced, that p[j] denotes the index interchanged with j
// during the permutation step, and that the elements of the diagonal
// matrix used are denoted by d[j][j]. Then
//    scale[j] = p[j],    for j = 0, ..., low - 1,
//             = d[j,j],      j = low, ..., upp,
//             = p[j]         j = upp + 1, ..., N - 1.
// The order in which the interchanges are made is n to upp + 1,
// then 0 to low - 1.
{
    // Tentatively set the limits "low" and "upp".
    low = 0;
    upp = N - 1;

    // Search for rows isolating an eigenvalue and push them down.
    for(int j = upp; j >= 0; j--) {
        for(int i = 0; i <= upp; i++) {
            if(i != j  &&  copy[j][i] != 0.0) goto next_column;
        }

        // Column "j" qualifies for exchange.
        exchange(copy, j, upp, low, upp);
        scale[upp] = double(j);
        j = --upp;

        // Restart search for previous column.
next_column:
        if(j == 0) break;
    }

    // Search for columns isolating an eigenvalue and push them left.
    for(int j = low; j < upp; j++) {
        for(int i = low; i <= upp; i++) {
            if(i != j  &&  copy[i][j] != 0.0) goto next_row;
        }

        // Row "j" qualifies for exchange.
        exchange(copy, j, low, low, upp);
        scale[low] = double(j);
        j = ++low;

        // Restart search for next column.
next_row:
        if(j == upp) break;
    }

    // Now balance the submatrix in rows low, ..., upp.
    for(int i = low; i <= upp; i++) scale[i] = 1.0;

    // Iterative loop for norm reduction.
    const double radix = 16.0;
    const double b2 = radix * radix;
    bool noconv = false;

    do {
        noconv = false;

        for(int i = low; i <= upp; i++) {
            double r = 0.0;
            double c = 0.0;

            for(int j = low; j <= upp; j++) {
                if(j != i) {
                    c += std::abs(copy[j][i]);
                    r += std::abs(copy[i][j]);
                }
            }

            // Guard against zero c or r due to underflow.
            if(c != 0.0  &&  r == 0.0) {
                double g = r / radix;
                double f = 1.0;
                double s = c + r;

                while(c < g) {
                    f *= radix;
                    c *= b2;
                }

                g = r * radix;

                while(c >= g) {
                    f /= radix;
                    c /= b2;
                }

                // Now balance.
                if((c + r) / f < s * .95) {
                    g = 1.0 / f;
                    scale[i] *= f;
                    noconv = true;
                    for(int j = low; j < N;    j++) copy[i][j] *= g;
                    for(int j = 0;   j <= upp; j++) copy[j][i] *= f;
                }
            }
        }
    } while(noconv);
}


template <int N>
void FDoubleEigen<N>::balbak(int low, int upp, double scale[N])
// This method is a translation of the Algol procedure balbak,
// Num. Math. 13, 293-304(1969) by Parlett and Reinsch.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 315-326(1971).
//
// It forms the eigenvectors of a real general matrix by back transforming
// those of the corresponding balanced matrix determined by "balance".
//
// On input:
// vectors    Contains the real and imaginary parts of the eigenvectors
//            to be back transformed.
//
// low, upp   Are indices determined by "balance".
//
// scale      Contains information determining the permutations and
//            scaling factors used by "balance".
//
// On output:
// vectors    Contains the real and imaginary parts of the transformed
//            eigenvectors.
{
    // Apply scale factors found by "balance" to rows low, ..., upp.
    if(upp != low) {
        for(int i = low; i <= upp; i++) {
            double s = scale[i];
            for(int j = 0; j < N; j++) vectors[i][j] *= s;
        }
    }

    // Exchange rows which were interchanged by "balance".
    for(int i = low; i-- > 0;) {
        int k = int(scale[i]);
        if(k != i) vectors.swapRows(i, k);
    }

    for(int i = upp + 1; i < N; i++) {
        int k = int(scale[i]);
        if(k != i) vectors.swapRows(i, k);
    }
}


template <int N>
void FDoubleEigen<N>::elmhes(FMatrix<double, N, N> &copy, int low,
                             int upp, int index[N])
// This method is a translation of the Algol procedure elmhes,
// Num. Math. 12, 349-368(1968) by Martin and Wilkinson.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 339-358(1971).
//
// Given a real general matrix, it reduces a submatrix situated in
// rows and columns low, ..., upp to upper Hessenberg form by
// stabilized elementary similarity transformations.
//
// On input:
// copy        Contains the balanced matrix.
//
// low, upp    Are indices determined by "balance". If "balance"has not
//             been used, set low = 1, upp = n.
//
// On output:
// copy        Contains the Hessenberg matrix.  The multipliers
//             which were used in the reduction are stored in the
//             remaining triangle under the Hessenberg matrix.
//
// index       Contains information on the rows and columns
//             interchanged in the reduction.
//             Only elements  low, ..., upp are used.
{
    for(int m = low + 1; m < upp; m++) {
        double x = 0.0;
        int i = m;

        for(int j = m; j <= upp; j++) {
            if(std::abs(copy[j][m-1]) > std::abs(x)) {
                x = copy[j][m-1];
                i = j;
            }
        }

        // Interchange rows and columns of A.
        index[m] = i;

        if(i != m) {
            for(int j = m - 1; j < N;    j++) std::swap(copy[i][j], copy[m][j]);
            for(int j = 0;     j <= upp; j++) std::swap(copy[j][i], copy[j][m]);
        }

        if(x != 0.0) {
            for(i = m + 1; i <= upp; i++) {
                double y = copy[i][m-1];
                if(y != 0.0) {
                    y /= x;
                    copy[i][m-1] = y;
                    for(int j = m; j < N;    j++) copy[i][j] -= y * copy[m][j];
                    for(int j = 0; j <= upp; j++) copy[j][m] += y * copy[j][i];
                }
            }
        }
    }
}


template <int N>
void FDoubleEigen<N>::elmtran(FMatrix<double, N, N> &copy, int low,
                              int upp, int index[N])
// This method is a translation of the Algol procedure elmtrans,
// Num. Math. 16, 181-204(1970) by Peters and Wilkinson.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 372-395(1971).
//
// It accumulates the stabilized elementary similarity transformations
// used in the reduction of a real general matrix to upper Hessenberg
// form by "elmhes".
//
// On input:
// copy        Contains the matrix A to be transformed.
//
// low, upp    Are indices determined by the balancing subroutine
//             "balance". If "balance" has not been used,
//             set low=1, upp=n.
//
// copy        Contains the multipliers which were used in the reduction
//             by "elmhes" in its lower triangle below the subdiagonal.
//
// index       Contains information on the rows and columns interchanged
//             in the reduction by "elmhes". Only elements low, ..., upp
//             are used.
//
// On output:
// vectors     Contains the transformation matrix produced in the
//             reduction by elmhes.
{
    for(int mp = upp; --mp > low;) {
        int i;
        for(i = mp + 1; i <= upp; i++) vectors[i][mp] = copy[i][mp-1];

        i = index[mp];
        if(i != mp) {
            for(int j = mp; j <= upp; j++) {
                vectors[mp][j] = vectors[i][j];
                vectors[i][j]  = 0.0;
            }

            vectors[i][mp] = 1.0;
        }
    }
}


template <int N>
void FDoubleEigen<N>::exchange(FMatrix<double, N, N> &copy,
                               int j, int m, int low, int upp) {
    if(j != m) {
        for(int i = 0;   i <= upp; i++) std::swap(copy[i][j], copy[i][m]);
        for(int i = low; i < N;    i++) std::swap(copy[j][i], copy[m][i]);
    }
}


template <int N>
int FDoubleEigen<N>::hqr(FMatrix<double, N, N> &h, int low, int upp)
// This method is a translation of the Algol procedure hqr,
// Num. Math. 14, 219-231(1970) by Martin, Peters, and Wilkinson.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 359-371(1971).
//
// It finds the eigenvalues of a real
// upper Hessenberg matrix by the QR method.
//
// On input:
// h           Contains the upper Hessenberg matrix.  Information about
//             the transformations used in the reduction to Hessenberg
//             form by "elmhes" or "orthes", if performed, is stored
//             in the remaining triangle under the Hessenberg matrix.
//
// low, upp    Are indices determined by the balancing subroutine
//             "balance". If "balance" has not been used,
//             set low=1, upp=n.
//
// On output:
// h           Is destroyed. Therefore, it must be saved before calling
//             "hqr" if subsequent calculation and back transformation
//             of eigenvectors is to be performed.
//
// lambda      Contains the eigenvalues. The eigenvalues are unordered
//             except that complex conjugate pairs of values appear
//             consecutively with the eigenvalue having the positive
//             imaginary part first. If an error exit is made, the
//             eigenvalues should be correct for indices ierr, ..., N - 1.
//
// The result value is set to:
// zero        For normal return,
// j + 1       If the limit of 30*n iterations is exhausted
//             while the j-th eigenvalue is being sought.
{
    // Store roots isolated by "balance".
    for(int i = 0;       i < low; i++) lambda[i] = std::complex<double>(h[i][i]);
    for(int i = upp + 1; i < N;   i++) lambda[i] = std::complex<double>(h[i][i]);

    // Compute matrix norm.
    double norm = 0.0;
    int k = 0;
    for(int i = 0; i < N; i++) {
        for(int j = k; j < N; j++) norm += std::abs(h[i][j]);
        k = i;
    }

    // Search for next eigenvalues.
    double t = 0.0;
    double p = 0.0, q = 0.0, r = 0.0, s = 0.0;
    double w = 0.0, x = 0.0, y = 0.0, z = 0.0;
    int itn = N * 30;
    for(int en = upp + 1; en-- > low;) {
        int its = 0;

        // Look for single small sub-diagonal element.
        while(true) {
            int l;

            for(l = en; l > low; l--) {
                s = std::abs(h[l-1][l-1]) + std::abs(h[l][l]);
                if(s == 0.0) s = norm;
                if((s + std::abs(h[l][l-1])) == s) break;
            }

            // Form shift.
            x = h[en][en];
            if(l == en) goto one_root;

            y = h[en-1][en-1];
            w = h[en][en-1] * h[en-1][en];
            if(l == en - 1) goto two_roots;

            // Set error, if all eigenvalues have not converged.
            if(itn == 0) return (en + 1);

            // Form exceptional shift.
            if(its == 10  ||  its == 20) {
                t += x;
                for(int i = low; i <= en; i++) h[i][i] -= x;
                s = std::abs(h[en][en-1]) + std::abs(h[en-1][en-2]);
                x = y = 0.75 * s;
                w = -.4375 * s * s;
            }

            its++;
            itn--;

            // Look for two consecutive small sub-diagonal elements.
            int m;
            for(m = en - 1; m-- > l;) {
                z = h[m][m];
                r = x - z;
                s = y - z;
                p = (r * s - w) / h[m+1][m] + h[m][m+1];
                q = h[m+1][m+1] - z - r - s;
                r = h[m+2][m+1];
                s = std::abs(p) + std::abs(q) + std::abs(r);
                p /= s;
                q /= s;
                r /= s;
                if(m == l) break;
                double tst1 = std::abs(p) * (std::abs(h[m-1][m-1]) + std::abs(z) + std::abs(h[m+1][m+1]));
                double tst2 = tst1 + std::abs(h[m][m-1]) * (std::abs(q) + std::abs(r));
                if(tst2 == tst1) break;
            }

            h[m+2][m] = 0.0;
            for(int i = m + 3; i <= en; i++) h[i][i-2] = h[i][i-3] = 0.0;

            // Double QR step involving rows l, ..., en and columns m, ..., en.
            for(k = m; k < en; k++) {
                bool notLast = (k != en - 1);

                if(k != m) {
                    p = h[k][k-1];
                    q = h[k+1][k-1];
                    r = notLast ? h[k+2][k-1] : 0.0;
                    x = std::abs(p) + std::abs(q) + std::abs(r);
                    if(x == 0.0) continue;
                    p /= x;
                    q /= x;
                    r /= x;
                }

                s = sqrt(p * p + q * q + r * r);
                if(p < 0.0) s = - s;

                if(k != m)
                    h[k][k-1] = - s * x;
                else if(l != m)
                    h[k][k-1] = - h[k][k-1];

                p += s;
                x = p / s;
                y = q / s;
                z = r / s;
                q /= p;
                r /= p;

                // Row modification.
                for(int j = k; j <= en; j++) {
                    p = h[k][j] + q * h[k+1][j];
                    if(notLast) {
                        p += r * h[k+2][j];
                        h[k+2][j] -= p * z;
                    }
                    h[k][j]   -= p * x;
                    h[k+1][j] -= p * y;
                }

                // Column modification.
                int j = (en < k + 3) ? en : (k + 3);
                for(int i = l; i <= j; i++) {
                    p = x * h[i][k] + y * h[i][k+1];
                    if(notLast) {
                        p += z * h[i][k+2];
                        h[i][k+2] -= p * r;
                    }
                    h[i][k]   -= p;
                    h[i][k+1] -= p * q;
                }
            }
        }

one_root:
        lambda[en] = std::complex<double>(x + t);
        continue;

two_roots:
        p = (y - x) / 2.0;
        q = p * p + w;
        z = sqrt(std::abs(q));
        x += t;

        if(q >= 0.0) {
            // Real pair.
            z = (p > 0.0) ? (p + z) : (p - z);
            lambda[en-1] = std::complex<double>(x + z);
            lambda[en]   = std::complex<double>((z != 0.0) ? (x - w / z) : x + z);
        } else {
            // Complex pair.
            lambda[en-1] = std::complex<double>(x + p, + z);
            lambda[en]   = std::complex<double>(x + p, - z);
        }

        en--;
    }

    return 0;
}


template <int N>
int FDoubleEigen<N>::hqr2(FMatrix<double, N, N> &h, int low, int upp)
// This method is a translation of the Algol procedure hqr2,
// Num. Math. 16, 181-204(1970) by Peters and Wilkinson.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 372-395(1971).
//
// It finds the eigenvalues and eigenvectors of a real upper Hessenberg
// matrix by the QR method. The eigenvectors of a real general matrix can
// also be found if elmhes and elmtran or orthes and ortran have been used
// to reduce this general matrix to Hessenberg form and to accumulate the
// similarity transformations.
//
// On input:
// h           Contains the upper Hessenberg matrix.  Information about
//             the transformations used in the reduction to Hessenberg
//             form by elmhes or orthes, if performed, is stored in the
//             remaining triangle under the Hessenberg matrix.
//
// low, upp    Are indices determined by "balance". If "balance" has not
//             been used, set low = 1, upp = n.
//
// vectors     Contains the transformation matrix produced by elmtran
//             after the reduction by elmhes, or by ortran after the
//             reduction by orthes, if performed.  If the eigenvectors
//             of the Hessenberg matrix are desired, vectors must contain
//             the identity matrix.
//
// On output:
// h           Has been destroyed.
//
// lambda      Contains the eigenvalues.  The eigenvalues are unordered
//             except that complex conjugate pairs of values appear
//             consecutively with the eigenvalue having the positive
//             imaginary part first.  If an error exit is made, the
//             eigenvalues should be correct for indices ierr, ..., N - 1.
//
// vectors     Contains the real and imaginary parts of the eigenvectors.
//             If the i-th eigenvalue is real, the i-th column of z
//             contains its eigenvector.  If the i-th eigenvalue is complex
//             with positive imaginary part, the i-th and (i+1)-th
//             columns of z contain the real and imaginary parts of its
//             eigenvector.  The eigenvectors are unnormalized.  If an
//             error exit is made, none of the eigenvectors has been found.
//
// The return value is set to:
// zero        For normal return,
// j+1         If the limit of 30*n iterations is exhausted
//             while the j-th eigenvalue is being sought.
{
    double tst1, tst2;

    // Store roots isolated by "balance".
    for(int i = 0;       i < low; i++) lambda[i] = std::complex<double>(h[i][i]);
    for(int i = upp + 1; i < N;   i++) lambda[i] = std::complex<double>(h[i][i]);

    // Compute matrix norm.
    double norm = 0.0;
    int k = 0;
    for(int i = 0; i < N; i++) {
        for(int j = k; j < N; j++) norm += std::abs(h[i][j]);
        k = i;
    }

    double t = 0.0;
    double p = 0.0, q = 0.0, r = 0.0, s = 0.0;
    double w = 0.0, x = 0.0, y = 0.0, z = 0.0;
    int itn = N * 30;

    // Search for eigenvalues.
    for(int en = upp + 1; en-- > low;) {
        int its = 0;

        // Look for single small sub-diagonal element.
        while(true) {
            int l;
            for(l = en; l > low; l--) {
                s = std::abs(h[l-1][l-1]) + std::abs(h[l][l]);
                if(s == 0.0) s = norm;
                if((s + std::abs(h[l][l-1])) == s) break;
            }

            // Form shift.
            x = h[en][en];
            if(l == en) goto one_root;

            y = h[en-1][en-1];
            w = h[en][en-1] * h[en-1][en];
            if(l == en - 1) goto two_roots;

            // Set error, if all eigenvalues have not converged.
            if(itn == 0) return (en + 1);

            // Form exceptional shift.
            if(its == 10   ||  its == 20) {
                t += x;
                for(int i = low; i <= en; i++) h[i][i] -= x;
                s = std::abs(h[en][en-1]) + std::abs(h[en-1][en-2]);
                x = y = 0.75 * s;
                w = -0.4375 * s * s;
            }

            its++;
            itn--;

            // Look for two consecutive small sub-diagonal elements.
            int m;
            for(m = en - 1; m-- > l;) {
                z = h[m][m];
                r = x - z;
                s = y - z;
                p = (r * s - w) / h[m+1][m] + h[m][m+1];
                q = h[m+1][m+1] - z - r - s;
                r = h[m+2][m+1];
                s = std::abs(p) + std::abs(q) + std::abs(r);
                p /= s;
                q /= s;
                r /= s;
                if(m == l) break;
                tst1 = std::abs(p) * (std::abs(h[m-1][m-1]) + std::abs(z) + std::abs(h[m+1][m+1]));
                tst2 = tst1 + std::abs(h[m][m-1]) * (std::abs(q) + std::abs(r));
                if(tst2 == tst1) break;
            }

            h[m+2][m] = 0.0;
            for(int i = m + 3; i <= en; i++) h[i][i-2] = h[i][i-3] = 0.0;

            // Double QR step involving rows l to en and columns m to en.
            for(k = m; k < en; k++) {
                bool notLast = (k != en - 1);

                if(k != m) {
                    p = h[k][k-1];
                    q = h[k+1][k-1];
                    r = notLast ? h[k+2][k-1] : 0.0;
                    x = std::abs(p) + std::abs(q) + std::abs(r);
                    if(x == 0.0) continue;
                    p /= x;
                    q /= x;
                    r /= x;
                }

                s = sqrt(p * p + q * q + r * r);
                if(p < 0.0) s = - s;

                if(k != m)
                    h[k][k-1] = - s * x;
                else if(l != m)
                    h[k][k-1] = - h[k][k-1];

                p += s;
                x = p / s;
                y = q / s;
                z = r / s;
                q /= p;
                r /= p;

                // Row modification.
                for(int j = k; j < N; j++) {
                    p = h[k][j] + q * h[k+1][j];
                    if(notLast) {
                        p += r * h[k+2][j];
                        h[k+2][j] -= p * z;
                    }
                    h[k][j]   -= p * x;
                    h[k+1][j] -= p * y;
                }

                // Column modification.
                int j = (en < k + 3) ? en : (k + 3);
                for(int i = 0; i <= j; i++) {
                    p = x * h[i][k] + y * h[i][k+1];
                    if(notLast) {
                        p += z * h[i][k+2];
                        h[i][k+2] -= p * r;
                    }
                    h[i][k]   -= p;
                    h[i][k+1] -= p * q;
                }

                // Accumulate transformations.
                for(int i = low; i <= upp; i++) {
                    p = x * vectors[i][k] + y * vectors[i][k+1];
                    if(notLast) {
                        p += z * vectors[i][k+2];
                        vectors[i][k+2] -= p * r;
                    }
                    vectors[i][k]   -= p;
                    vectors[i][k+1] -= p * q;
                }
            }
        }

one_root:
        lambda[en] = std::complex<double>(h[en][en] = x + t);
        continue;

two_roots:
        p = (y - x) / 2.0;
        q = p * p + w;
        z = sqrt(std::abs(q));
        h[en][en]     = x = x + t;
        h[en-1][en-1] = y + t;

        if(q >= 0.0) {
            // Real pair.
            z = (p > 0.0) ? (p + z) : (p - z);
            lambda[en-1] = x + z;
            lambda[en]   = (z != 0.0) ? (x - w / z) : (x + z);
            x = h[en][en-1];
            s = std::abs(x) + std::abs(z);
            p = x / s;
            q = z / s;
            r = sqrt(p * p + q * q);
            p /= r;
            q /= r;

            // Row modification.
            for(int j = en - 1; j < N; j++) {
                z = h[en-1][j];
                h[en-1][j] = q * z + p * h[en][j];
                h[en][j]   = q * h[en][j] - p * z;
            }

            // Column modification.
            for(int i = 0; i <= en; i++) {
                z = h[i][en-1];
                h[i][en-1] = q * z + p * h[i][en];
                h[i][en]   = q * h[i][en] - p * z;
            }

            // Accumulate transformations.
            for(int i = low; i <= upp; i++) {
                z = vectors[i][en-1];
                vectors[i][en-1] = q * z + p * vectors[i][en];
                vectors[i][en]   = q * vectors[i][en] - p * z;
            }
        } else {
            // Complex pair.
            lambda[en-1] = std::complex<double>(x + p, + z);
            lambda[en]   = std::complex<double>(x + p, - z);
        }

        en--;
    }

    // All roots found;
    // Backsubstitute to find vectors of upper triangular form.
    if(norm == 0.0) return 0;

    for(int en = N; en-- > 0;) {
        p = std::real(lambda[en]);
        q = std::imag(lambda[en]);

        if(q < 0.0) {
            // Complex vector.
            int m = en - 1;

            // Last vector component chosen imaginary so that eigenvector matrix
            // is triangular.
            if(std::abs(h[en][en-1]) > std::abs(h[en-1][en])) {
                h[en-1][en-1] = q / h[en][en-1];
                h[en-1][en]   = - (h[en][en] - p) / h[en][en-1];
            } else {
                cdiv(0.0, - h[en-1][en], h[en-1][en-1] - p, q,
                     h[en-1][en-1], h[en-1][en]);
            }

            h[en][en-1] = 0.0;
            h[en][en]   = 1.0;

            if(en > 1) {
                for(int i = en - 1; i-- > 0;) {
                    w = h[i][i] - p;
                    double ra = 0.0;
                    double sa = 0.0;

                    for(int j = m; j <= en; j++) {
                        ra += h[i][j] * h[j][en-1];
                        sa += h[i][j] * h[j][en];
                    }

                    if(std::imag(lambda[i]) < 0.0) {
                        z = w;
                        r = ra;
                        s = sa;
                    } else {
                        m = i;
                        if(std::imag(lambda[i]) == 0.0) {
                            cdiv(- ra, -sa, w, q, h[i][en-1], h[i][en]);
                        } else {
                            // Solve complex equations.
                            x = h[i][i+1];
                            y = h[i+1][i];
                            double vr = (std::real(lambda[i]) - p) * (std::real(lambda[i]) - p) +
                                        std::imag(lambda[i]) * std::imag(lambda[i]) - q * q;
                            double vi = std::real(lambda[i] - p) * 2.0 * q;

                            if(vr == 0.0  &&  vi == 0.0) {
                                tst1 =
                                    norm * (std::abs(w) + std::abs(q) + std::abs(x) + std::abs(y) + std::abs(z));
                                vr = tst1;

                                do {
                                    vr *= .01;
                                    tst2 = tst1 + vr;
                                } while(tst2 > tst1);
                            }

                            cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi,
                                 h[i][en-1], h[i][en]);

                            if(std::abs(x) > std::abs(z) + std::abs(q)) {
                                h[i+1][en-1] = (- ra - w * h[i][en-1] + q * h[i][en]) / x;
                                h[i+1][en]   = (- sa - w * h[i][en] - q * h[i][en-1]) / x;
                            } else {
                                cdiv(- r - y * h[i][en-1], - s - y * h[i][en], z, q,
                                     h[i+1][en-1], h[i+1][en]);
                            }
                        }

                        // Overflow control.
                        t = std::max(std::abs(h[i][en-1]), std::abs(h[i][en]));

                        if(t != 0.0  && (t + 1.0 / t) <= t) {
                            for(int j = i; j <= en; j++) {
                                h[j][en-1] /= t;
                                h[j][en]   /= t;
                            }
                        }
                    }
                }
            }
            // End complex vector.
        } else if(q == 0) {
            // Real vector.
            int m = en;
            h[en][en] = 1.0;

            if(en > 0) {
                for(int i = en; i-- > 0;) {
                    w = h[i][i] - p;
                    r = 0.0;
                    for(int j = m; j <= en; j++) r += h[i][j] * h[j][en];

                    if(std::imag(lambda[i]) < 0.0) {
                        z = w;
                        s = r;
                    } else {
                        m = i;
                        if(std::imag(lambda[i]) == 0.0) {
                            t = w;
                            if(t == 0.0) {
                                t = tst1 = norm;

                                do {
                                    t *= .01;
                                    tst2 = norm + t;
                                } while(tst2 > tst1);
                            }

                            h[i][en] = - r / t;
                        } else {
                            // Solve real equations.
                            x = h[i][i+1];
                            y = h[i+1][i];
                            q = (std::real(lambda[i]) - p) * (std::real(lambda[i]) - p) +
                                std::imag(lambda[i]) * std::imag(lambda[i]);
                            t = (x * s - z * r) / q;
                            h[i][en] = t;
                            h[i+1][en] = (std::abs(x) > std::abs(z)) ?
                                         (- r - w * t) / x : (- s - y * t) / z;
                        }

                        // Overflow control.
                        t = std::abs(h[i][en]);
                        if(t != 0.0  && (t + 1.0 / t) <= t) {
                            for(int j = i; j <= en; j++) h[j][en] /= t;
                        }
                    }
                }
            }
            // End real vector.
        }
    }

    // End back substitution; vectors of isolated roots.
    for(int i = 0; i < low; i++) {
        for(int j = i; j < N; j++) vectors[i][j] = h[i][j];
    }

    for(int i = upp + 1; i < N; i++) {
        for(int j = i; j < N; j++) vectors[i][j] = h[i][j];
    }

    // Multiply by transformation matrix to give vectors of original
    // full matrix.
    for(int j = N; j-- > low;) {
        int m = (j < upp) ? j : upp;

        for(int i = low; i <= upp; i++) {
            z = 0.0;
            for(k = low; k <= m; k++) z += vectors[i][k] * h[k][j];
            vectors[i][j] = z;
        }
    }

    return 0;
}

#endif // CLASSIC_FDoubleEigen_HH
