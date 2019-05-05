// ------------------------------------------------------------------------
// $RCSfile: ComplexEigen.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class ComplexEigen:
//   Eigenvector methods for class: Matrix<complex<double> >.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/ComplexEigen.h"
#include "Algebra/Matrix.h"
#include "Utilities/EigenvalueError.h"
#include "Utilities/SizeError.h"
#include <cmath>

using std::abs;
using std::conj;
using std::imag;
using std::norm;
using std::real;
using std::sqrt;
using std::swap;

// Helper routine.
// ------------------------------------------------------------------------

namespace {
    inline double abssum(complex<double> a) {
        return std::abs(real(a)) + std::abs(imag(a));
    }
}


// Class ComplexEigen
// --------------------------------------------------------------------------

ComplexEigen::ComplexEigen():
    lambda(), vectors()
{}


ComplexEigen::ComplexEigen(const ComplexEigen &rhs):
    lambda(rhs.lambda), vectors(rhs.vectors)
{}


ComplexEigen::ComplexEigen(const Matrix<complex<double> > &M, bool vec):
    lambda(M.nrows()), vectors(M.nrows(), M.nrows(), 1.0) {
    int n = M.nrows();
    int m = M.ncols();
    if(n != m) {
        throw SizeError("ComplexEigen::ComplexEigen()", "Matrix is not square.");
    }
    Matrix<complex<double> > copy(M);

    int low, upp;
    Array1D<double> scale(n);
    balance(copy, low, upp, scale);

    Array1D<complex<double> > ort(n);
    orthes(copy, low, upp, ort);

    if(vec) {
        if(hqr2(copy, low, upp, ort) != 0) {
            throw EigenvalueError("ComplexEigen::ComplexEigen()",
                                  "Unable to find all eigenvalues.");
        } else {
            balbak(low, upp, scale);
        }
    } else {
        if(hqr(copy, low, upp) != 0) {
            throw EigenvalueError("ComplexEigen::ComplexEigen()",
                                  "Unable to find all eigenvalues.");
        }
    }
}


ComplexEigen::~ComplexEigen()
{}


const Vector<complex<double> > &ComplexEigen::eigenValues() const {
    return lambda;
}


const Matrix<complex<double> > &ComplexEigen::eigenVectors() const {
    return vectors;
}


void ComplexEigen::balance(Matrix<complex<double> > &copy, int &low,
                           int &upp, Array1D<double> &scale)
// This subroutine is a translation of the Algol procedure "cbalance",
// a complex version of "balance",
// Num. Math. 13, 293-304(1969) by Parlett and Reinsch.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 315-326(1971).
//
// This subroutine balances a complex matrix and isolates eigenvalues
// whenever possible.
//
// On input:
// copy       Contains the complex<double> matrix to be balanced.
//
// On output:
// copy       Contains the balanced matrix.
//
// low, upp   Are two integers such that A[i][j] is zero if:
//            (1) i is greater than j and
//            (2) j = 0, ..., low - 1 or i = upp + 1, ..., n - 1.
//
// scale      Contains information determining the permutations and
//            scaling factors used.
//
// Suppose that the principal submatrix in rows low through upp has
// been balanced, that p[j] denotes the index interchanged with j
// during the permutation step, and that the elements of the diagonal
// matrix used are denoted by d[i,j]. Then
//    scale[j] = p[j],    for j = 0, ..., low - 1
//             = d[j,j]       j = low, ..., upp
//             = p[j]         j = upp + 1, ..., n - 1.
// The order in which the interchanges are made is n - 1 to upp + 1,
// then 0 to low - 1.
//
// Arithmetic is complex<double> throughout.
//
// Translated to FORTRAN by Burton S. Garbov, ANL
// Adapted March 1997 by F. Christoph Iselin, CERN, SL/AP
{
    int n = copy.nrows();
    low = 0;
    upp = n - 1;

    // Search for rows isolating an eigenvalue and push them down.
    for(int j = upp; ; j--) {
        for(int i = 0; i <= upp; i++) {
            if(i != j  &&  copy[j][i] != 0.0) goto next_row;
        }

        // Column "j" qualifies for exchange.
        exchange(copy, j, upp, low, upp);
        scale[upp] = double(j);
        j = --upp;

        // Restart search for next column.
next_row:
        if(j == 0) break;
    }

    // Search for columns isolating an eigenvalue and push them left.
    for(int j = low; ; j++) {
        for(int i = low; i <= upp; i++) {
            if(i != j  &&  copy[i][j] != 0.0) goto next_column;
        }

        // Column "j" qualifies for exchange.
        exchange(copy, j, low, low, upp);
        scale[low] = double(j);
        j = ++low;

        // Restart search for next column.
        /*
        Die Abbruchbedingung in der Methode balance auf Zeile 176 der Datei /scratch2/amas/l_felsimsvn/src/opal/classic/5.0/src/Algebra/ComplexEigen.cpp
        sollte von
         if(j = upp) break;
        auf
         if(j >= upp) break;
        geändert werden.

        In der bisherige Version wird die Abbruchbedingung in seltenen Fällen überhüpft, was zu einer unendlichen Schleife und einem unkontrollierten Programmabsturz führt.
        Dieser Fall tritt auf, wenn die analysierte Matrix schon diagonal ist, z.B. bei folgendem einfachen Beispiel:

         Matrix<double> m(2,2,0.0);

         m[0][0]=1;
         m[1][1]=3;

         DoubleEigen de(m,true);

         Dies ist zwar ein bisschen ein akademischer Fall, war aber ausgerechnet mein erstes Testbeispiel.
         Die vorgeschlagene Änderung ist performancemässig gleichwertig und sollte keine negativen Seiteneffekte haben.

        */
next_column:
        if(j >= upp) break;
    }

    // Initialize scale factors.
    for(int i = low; i <= upp; i++) scale[i] = 1.0;

    // Now balance the submatrix in rows low to upp.
    const double radix = 16.;
    const double b2 = radix * radix;
    bool noconv;

    // Iterative loop for norm reduction.
    do {
        noconv = false;

        for(int i = low; i <= upp; i++) {
            double c = 0.0;
            double r = 0.0;

            for(int j = low; j <= upp; j++) {
                if(j != i) {
                    c = c + abssum(copy[j][i]);
                    r = r + abssum(copy[i][j]);
                }
            }

            // Guard against zero c or r due to underflow/
            if(c != 0.0  &&  r != 0.0) {
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

                // Now balance/
                if((c + r) / f < s * .95) {
                    g = 1.0 / f;
                    scale[i] *= f;
                    noconv = true;
                    for(int j = low; j <= n; j++) copy[i][j] *= g;
                    for(int j = 1; j <= upp; j++) copy[j][i] *= f;
                }
            }
        }
    } while(noconv);

    return;
}


void ComplexEigen::balbak(int low, int upp,
                          const Array1D<double> &scale)
// This subroutine is a translation of the Algol procedure cbabk2,
// a complex<double> version of balbak,
// Num. Math. 13, 293-304(1969) by Parlett and Reinsch.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 315-326(1971).
//
// It forms the eigenvectors of a complex general matrix by back
// transforming those of the corresponding balanced matrix determined
// by "balance".
//
// On input:
// vectors   Contains the eigenvectors to be back transformed.
//
// low, upp  Are integers determined by "balance".
//
// scale     Contains information determining the permutations
//           and scaling factors used by  "balance".
//
// On output:
// vectors   Contains the the transformed eigenvectors.
//
// Translated to FORTRAN by Burton S. Garbov, ANL
// Adapted March 1997 by F. Christoph Iselin, CERN, SL/AP
{
    int n = vectors.nrows();

    // Apply scale factors found by "balance" to rows low, ..., upp.
    if(upp != low) {
        for(int i = low; i <= upp; i++) {
            double s = scale[i];
            for(int j = 0; j < n; j++) vectors[i][j] *= s;
        }
    }

    // Exchange rows which were interchanged by "balance".
    // for i = low-1 step -1 until 0, upp+1 step 1 until n-1 do
    for(int i = low; i-- > 0;) {
        int k = int(scale[i]);
        if(k != i) vectors.swapRows(k, i);
    }

    for(int i = upp + 1; i < n; i++) {
        int k = int(scale[i]);
        if(k != i) vectors.swapRows(k, i);
    }

    return;
}


void ComplexEigen::exchange
(Matrix<complex<double> > &copy, int j, int m, int low, int upp) {
    int n = copy.nrows();

    if(j != m) {
        for(int i = 0; i <= upp; i++) swap(copy[i][j], copy[i][m]);
        for(int i = low; i <= n; i++) swap(copy[j][i], copy[m][i]);
    }
}


int ComplexEigen::hqr(Matrix<complex<double> > &h, int low, int upp)
// This subroutine is a translation of a unitary analogue of the
// Algol procedure  comlr, Num. Math. 12, 369-376(1968) by Martin
// and Wilkinson.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 396-403(1971).
// the unitary analogue substitutes the QR algorithm of Francis
// (Comp. Jour. 4, 332-345(1962)) for the LR algorithm.
//
// This subroutine finds the eigenvalues of a complex
// upper Hessenberg matrix by the QR method.
//
// On input:
// h         Contains the complex<double> upper Hessenberg matrix.
//           Its lower triangles below the subdiagonal contain
//           information about the unitary transformations used in
//           the reduction by  orthes, if performed.
//
// low, upp  Are integers determined by the balancing subroutine "balance".
//           if  "balance"  has not been used, set low=1, upp=n.
//
// On output:
// h         The upper Hessenberg portion has been destroyed.
//           Therefore, it must be saved before calling hqr
//           if subsequent calculation of eigenvectors is to be performed.
//
// lambda    Contains the eigenvalues.  If an error exit is made,
//           the eigenvalues should be correct for indices ierr,...,n.
//
// The return value is
// zero      For normal return,
// j+1       If the limit of 30n iterations is exhausted
//           while the j-th eigenvalue is being sought.
//
// Translated to FORTRAN by Burton S. Garbov, ANL
// Adapted March 1997 by F. Christoph Iselin, CERN, SL/AP
{
    int n = h.nrows();
    complex<double> s, x, y, z;

    // Create real subdiagonal elements.
    for(int i = low + 1; i <= upp; i++) {
        if(imag(h[i][i-1]) != 0.0) {
            double norm = std::abs(h[i][i-1]);
            y = h[i][i-1] / norm;
            h[i][i-1] = complex<double>(norm, 0.0);
            for(int j = i; j <= upp; j++) h[i][j] = conj(y) * h[i][j];
            int ll = (upp <= i) ? upp : (i + 1);
            for(int j = low; j <= ll; j++) h[j][i] = y * h[j][i];
        }
    }

    // Store roots isolated by "balance".
    for(int i = 0; i < n; i++) {
        if(i < low  ||  i < upp) lambda[i] = h[i][i];
    }

    // Search for eigenvalues,
    complex<double> t = 0.0;
    int itn = n * 30;

    for(int en = upp + 1; en-- > low;) {
        int its = 0;

        // Look for single small sub-diagonal element.
        while(true) {
            int l;
            for(l = en; l > low; l--) {
                double tst1, tst2;
                tst1 = abssum(h[l-1][l-1]) + abssum(h[l][l]);
                tst2 = tst1 + std::abs(real(h[l][l-1]));
                if(tst2 == tst1) break;
            }

            // Form shift.
            if(l == en) break;

            // Set error -- all eigenvalues have not converged after 30n iterations.
            if(itn == 0) return (en + 1);

            if(its != 10  &&  its != 20) {
                s = h[en][en];
                x = h[en-1][en] * real(h[en][en-1]);

                if(x != 0.0) {
                    y = (h[en-1][en-1] - s) / 2.0;
                    z = sqrt(y * y + x);

                    if(real(y) * real(z) + imag(y) * imag(z) < 0.0)
                        z = - z;

                    x /= (y + z);
                    s -= x;
                }
            } else {
                // Form exceptional shift.
                s = std::abs(real(h[en][en-1])) + std::abs(real(h[en-1][en-2]));
            }

            for(int i = low; i <= en; i++) h[i][i] -= s;
            t += s;
            its++;
            itn--;

            // Reduce to triangle (rows).
            for(int i = l + 1; i <= en; i++) {
                double sr = real(h[i][i-1]);
                double norm = hypot(std::abs(h[i-1][i-1]), sr);
                lambda[i-1] = x = h[i-1][i-1] / norm;
                h[i-1][i-1] = norm;
                double fi = sr / norm;
                h[i][i-1] = complex<double>(0.0, fi);

                for(int j = i; j <= en; j++) {
                    y = h[i-1][j];
                    z = h[i][j];
                    h[i-1][j] = conj(x) * y + fi * z;
                    h[i][j]   = x * z       - fi * y;
                }
            }

            double si = imag(h[en][en]);

            if(si != 0.0) {
                double norm = std::abs(h[en][en]);
                s = h[en][en] / norm;
                h[en][en] = norm;
            }

            // Inverse operation (columns).
            for(int j = l + 1; j <= en; j++) {
                double fi = imag(h[j][j-1]);
                x = lambda[j-1];

                for(int i = l; i < j; i++) {
                    y = h[i][j-1];
                    z = h[i][j];
                    h[i][j-1] = x * y       + fi * z;
                    h[i][j]   = conj(x) * z - fi * y;
                }

                double yr   = real(h[j][j-1]);
                z           = h[j][j];
                h[j][j-1] = real(x) * yr + fi * real(z);
                h[j][j]   = conj(x) * z  - fi * yr;
            }

            if(si != 0.0) {
                for(int i = l; i <= en; i++) h[i][en] *= s;
            }
        }

        // A root found.
        lambda[en] = h[en][en] + t;
    }

    // All eigenvalues have been found.
    return 0;
}


int ComplexEigen::hqr2(Matrix<complex<double> > &h, int low,
                       int upp, Array1D<complex<double> > &ort)
// This subroutine is a translation of a unitary analogue of the
// Algol procedure  comlr2, Num. Math. 16, 181-204(1970) by Peters
// and Wilkinson.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 372-395(1971).
// the unitary analogue substitutes the QR algorithm of Francis
// (Comp. Jour. 4, 332-345(1962)) for the LR algorithm.
//
// This subroutine finds the eigenvalues and eigenvectors
// of a complex<double> upper Hessenberg matrix by the qr
// method.  The eigenvectors of a complex<double> general matrix
// can also be found if  orthes  has been used to reduce
// this general matrix to Hessenberg form.
//
// On input:
// h         Contains the complex<double> upper Hessenberg matrix.
//           Its lower triangle below the subdiagonal contains further
//           information about the transformations which were used in the
//           reduction by  orthes, if performed.  If the eigenvectors of
//           the Hessenberg matrix are desired, these elements may be
//           arbitrary.
//
// low, upp  Are integers determined by the balancing subroutine "balance".
//           If  "balance"  has not been used, set low=1, upp=n.
//
// ort       Contains information about the unitary transformations
//           used in the reduction by  orthes, if performed.
//           Only elements low through upp are used.  If the eigenvectors
//           of the Hessenberg matrix are desired, set ort(j) to zero for
//           these elements.
//
// On output:
// ort       And the upper Hessenberg portion of  h  have been destroyed.
//
// lambda    Contains the eigenvalues.  if an error exit is made,
//           the eigenvalues should be correct for indices ierr,...,n.
//
// vectors   Contains the eigenvectors.  The eigenvectors are unnormalized.
//           If an error exit is made, none of the eigenvectors has been
//           found.
//
// The return value is
// zero      For normal return,
// j+1       If the limit of 30n iterations is exhausted
//           while the j-th eigenvalue is being sought.
//
// Translated to FORTRAN by Burton S. Garbov, ANL
// Adapted March 1997 by F. Christoph Iselin, CERN, SL/AP
{
    ssize_t n = h.nrows();
    complex<double> s, x, y, z;

    // Form the matrix of accumulated transformations from the information
    // left by "orthes".
    for(ssize_t i = upp - 1; i > low; i--) {
        if(ort[i] != 0.0  &&  h[i][i-1] != 0.0) {
            // Norm below is negative of h formed in orthes
            double norm = real(h[i][i-1]) * real(ort[i]) +
                          imag(h[i][i-1]) * imag(ort[i]);

            for(ssize_t k = i + 1; k <= upp; k++) ort[k] = h[k][i-1];

            for(ssize_t j = i; j <= upp; j++) {
                s = 0.0;
                for(ssize_t k = i; k <= upp; k++) s += conj(ort[k]) * vectors[k][j];
                s /= norm;
                for(ssize_t k = i; k <= upp; k++) vectors[k][j] += s * ort[k];
            }
        }
    }

    // Create real subdiagonal elements.
    for( ssize_t i = low + 1; i <= upp; i++) {
        if(imag(h[i][i-1]) != 0.0) {
            double norm = std::abs(h[i][i-1]);
            y = h[i][i-1] / norm;
            h[i][i-1] = norm;
            for(ssize_t j = i; j < n; j++) h[i][j] = conj(y) * h[i][j];
            int ll = (upp <= i) ? upp : (i + 1);
            for(ssize_t j = 0; j <= ll; j++) h[j][i] = y * h[j][i];
            for(ssize_t j = low; j <= upp; j++) vectors[j][i] = y * vectors[j][i];
        }
    }

    // Store roots isolated by "balance".
    for(ssize_t i = 0; i < n; i++) {
        if(i < low  ||  i > upp) lambda[i] = h[i][i];
    }

    complex<double> t = 0.0;
    ssize_t itn = n * 30;

    // Search for eigenvalues.
    for(ssize_t en = upp + 1; en-- > low;) {
        ssize_t its = 0;

        ssize_t l = en;
        // Look for singlengle small sub-diagonal element.
        while(true) {
            for(; l > low; l--) {
                double tst1, tst2;
                tst1 = abssum(h[l-1][l-1]) + abssum(h[l][l]);
                tst2 = tst1 + std::abs(real(h[l][l-1]));
                if(tst2 == tst1) break;
            }

            if(l == en) break;

            // Set error -- all eigenvalues have not converged after 30n iterations.
            if(itn == 0) return (en + 1);

            if(its != 10  &&  its != 20) {
                // Form shift.
                s = h[en][en];
                x = h[en-1][en] * real(h[en][en-1]);

                if(x != 0.0) {
                    y = (h[en-1][en-1] - s) / 2.0;
                    z = sqrt(y * y + x);

                    if(real(y) * real(z) + imag(y) * imag(z) < 0.0)
                        z = - z;

                    x /= (y + z);
                    s -= x;
                }
            } else {
                // Form exceptional shift.
                s = std::abs(real(h[en][en-1])) + std::abs(real(h[en-1][en-2]));
            }

            for(ssize_t i = low; i <= en; i++) h[i][i] -= s;
            t += s;
            its++;
            itn--;

            // Reduce to triangle (rows).
            for(ssize_t i = l + 1; i <= en; i++) {
                double sr = real(h[i][i-1]);
                double norm = hypot(std::abs(h[i-1][i-1]), sr);
                lambda[i-1] = x = h[i-1][i-1] / norm;
                h[i-1][i-1] = norm;
                double fi = sr / norm;
                h[i][i-1] = complex<double>(0.0, fi);

                for(ssize_t j = i; j < n; j++) {
                    y = h[i-1][j];
                    z = h[i][j];
                    h[i-1][j] = conj(x) * y + fi * z;
                    h[i][j]   = x       * z - fi * y;
                }
            }

            double si = imag(h[en][en]);

            if(si != 0.0) {
                double norm = std::abs(h[en][en]);
                s = h[en][en] / norm;
                h[en][en] = norm;
                for(ssize_t j = en + 1; j < n; j++) h[en][j] *= conj(s);
            }

            // Inverse operation (columns).
            for(ssize_t j = l + 1; j <= en; j++) {
                x = lambda[j-1];
                double fi = imag(h[j][j-1]);

                for(ssize_t i = 0; i < j; i++) {
                    y = h[i][j-1];
                    z = h[i][j];
                    h[i][j-1] = x       * y + fi * z;
                    h[i][j]   = conj(x) * z - fi * y;
                }

                double yr   = real(h[j][j-1]);
                z           = h[j][j];
                h[j][j-1] = real(x) * yr + fi * real(z);
                h[j][j]   = conj(x) * z  - fi * yr;

                for(ssize_t i = low; i <= upp; i++) {
                    y = vectors[i][j-1];
                    z = vectors[i][j];
                    vectors[i][j-1] = x       * y + fi * z;
                    vectors[i][j]   = conj(x) * z - fi * y;
                }
            }

            if(si != 0.0) {
                for(ssize_t i = 0; i <= en; i++) h[i][en] *= s;
                for(ssize_t i = low; i <= upp; i++) vectors[i][en] *= s;
            }
        }

        // A root found.
        h[en][en] += t;
        lambda[en] = h[en][en];
    }

    // All roots found.
    // Backsubstitute to find vectors of upper triangular form.
    double norm = 0.0;

    for(ssize_t i = 0; i < n; i++) {
        for(ssize_t j = i; j < n; j++) {
            double temp = abssum(h[i][j]);
            if(temp > norm) norm = temp;
        }
    }

    if(n == 1  ||  norm == 0.0) return 0;

    for(ssize_t en = n - 1; en > 0; en--) {
        x = lambda[en];
        h[en][en] = 1.0;

        for(ssize_t i = en; i-- > 0;) {
            z = 0.0;
            for(ssize_t j = i + 1; j <= en; j++) z += h[i][j] * h[j][en];
            y = x - lambda[i];

            if(y == 0.0) {
                double tst1 = norm;
                double tst2;

                do {
                    tst1 *= .01;
                    tst2 = norm + tst1;
                } while(tst2 > norm);

                y = tst1;
            }

            h[i][en] = z / y;

            // Overflow control.
            double temp = abssum(h[i][en]);

            if(temp != 0.0) {
                double tst1 = temp;
                double tst2 = tst1 + 1.0 / tst1;

                if(tst2 <= tst1) {
                    for(ssize_t j = i; j <= en; j++) h[j][en] /= temp;
                }
            }
        }
    }

    // End backsubstitution; vectors of isolated roots.
    for(ssize_t i = 0; i < n; i++) {
        if(i < low  ||  i > upp) {
            for(ssize_t j = i; j <= n; j++) vectors[i][j] = h[i][j];
        }
    }

    // Multiply by transformation matrix to give vectors of original full matrix.
    for(ssize_t j = n; j-- > low;) {
        ssize_t m = (j < upp) ? j : upp;

        for(ssize_t i = low; i <= upp; i++) {
            z = 0.0;
            for(ssize_t k = low; k <= m; k++) z += vectors[i][k] * h[k][j];
            vectors[i][j] = z;
        }
    }

    return 0;
}


void ComplexEigen::orthes(Matrix<complex<double> > &copy, int low,
                          int upp, Array1D<complex<double> > &ort)
// This subroutine is a translation of a complex<double> analogue of
// the Algol procedure orthes, Num. Math. 12, 349-368(1968)
// by Martin and Wilkinson.
// Handbook for Auto. Comp., Vol.ii-Linear Algebra, 339-358(1971).
//
// Given a complex<double> general matrix, this subroutine
// reduces a submatrix situated in rows and columns
// low through upp to upper Hessenberg form by
// unitary similarity transformations.
//
// On input:
// copy      Contains the complex<double> input matrix.
//
// low, upp  Are integers determined by the balancing subroutine "balance".
//           if  "balance"  has not been used, set low=1, upp=n.
//
// On output:
// copy      Contains the Hessenberg matrix.  Information about the
//           unitary transformations used in the reduction is stored in
//           the remaining triangles under the Hessenberg matrix.
//
// ort       Contains further information about the transformations.
//           Only elements low through upp are used.
//
// Translated to FORTRAN by Burton S. Garbov, ANL
// Adapted March 1997 by F. Christoph Iselin, CERN, SL/AP
{
    int n = copy.nrows();

    for(int m = low + 1; m < upp; m++) {
        double h = 0.0;
        ort[m] = 0.0;
        double scale = 0.0;

        // Scale column (Algol tol then not needed).
        for(int i = m; i <= upp; i++) scale += abssum(copy[i][m-1]);

        if(scale != 0.0) {
            for(int i = upp + 1; i-- > m;) {
                ort[i] = copy[i][m-1] / scale;
                h += norm(ort[i]);
            }

            double g = sqrt(h);
            double f = std::abs(ort[m]);

            if(f != 0.0) {
                h += f * g;
                g /= f;
                ort[m] *= (g + 1.0);
            } else {
                ort[m] = g;
                copy[m][m-1] = scale;
            }

            // Form (I - (u*ut)/h) * A.
            for(int j = m; j < n; j++) {
                complex<double> f = 0.0;
                for(int i = upp + 1; i-- > m;) f += conj(ort[i]) * copy[i][j];
                f /= h;
                for(int i = m; i <= upp; i++) copy[i][j] -= f * ort[i];
            }

            // Form (I - (u*ut)/h) * A * (I - (u*ut)/h).
            for(int i = 0; i <= upp; i++) {
                complex<double> f = 0.0;
                for(int j = upp + 1; j-- > m;) f += ort[j] * copy[i][j];
                f /= h;
                for(int j = m; j <= upp; j++) copy[i][j] -= f * conj(ort[j]);
            }

            ort[m] *= scale;
            copy[m][m-1] = - g * copy[m][m-1];
        }
    }
}