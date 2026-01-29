#ifndef CLASSIC_LUMatrix_HH
#define CLASSIC_LUMatrix_HH

// ------------------------------------------------------------------------
// $RCSfile: LUMatrix.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: LUMatrix
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"
#include "Algebra/Matrix.h"
#include "Utilities/SingularMatrixError.h"
#include "Utilities/SizeError.h"
#include <cmath>


// Template class LUMatrix<T>
// ------------------------------------------------------------------------
/// Triangle decomposition of a matrix.
//  A templated representation of a LU-decomposition.
//  Constructed from a Matrix<T>, this class remembers the triangular
//  decomposition and can return the back-substitution or the inverse.
//  When LUMatrix is used on a complex matrix, then the statement
//  [pre]
//     include std::abs;
//  [/pre]
//  is required in the calling file, so as to inject the [tt]abs()[/tt]
//  function to the global namespace.

template <class T>
class LUMatrix {

public:

    /// Constructor.
    //  Construct triangular decomposition of [b]m[/b].
    //  Throw [b]SingularMatrixError[/b] if the [b]m[/b] is singular.
    explicit LUMatrix(const Matrix<T> &m);

    LUMatrix();
    LUMatrix(const LUMatrix<T> &);
    ~LUMatrix();
    LUMatrix<T> &operator=(const LUMatrix<T> &);

    /// Back substitution.
    //  Perform back substitution on the vector [b]B[/b].
    //  The new [b]B[/b] is the old [b]B[/b], pre-multiplied by the inverse.
    void backSubstitute(Vector<T> &B) const;

    /// Back substitution.
    //  Perform back substitution on the matrix [b]M[/b].
    //  The new [b]M[/b] is the old [b]M[/b], pre-multiplied by the inverse.
    void backSubstitute(Matrix<T> &M) const;

    /// Get inverse.
    //  Construct and return the inverse.
    Matrix<T> inverse() const;

private:

    // Back substitution on iterator.
    template <class Iterator> void backSubstitute(Iterator) const;

    // The LU-decomposed matrix.
    Matrix<T> decomp;

    // The [b]index[/b] array used for LU decomposition and back substitution.
    Array1D<int> index;
};


// Template function implementation.
// ------------------------------------------------------------------------

template <class T>
LUMatrix<T>::LUMatrix():
    decomp(), index()
{}


template <class T>
LUMatrix<T>::LUMatrix(const LUMatrix<T> &rhs):
    decomp(rhs.decomp), index(rhs.index)
{}


template <class T>
LUMatrix<T>::LUMatrix(const Matrix<T> &rhs):
    decomp(rhs), index(rhs.nrows()) {
    if(rhs.nrows() != rhs.ncols()) {
        throw SizeError("LUMatrix::LUMatrix()", "Matrix is not square.");
    }

    // Builds the L-U decomp of a square matrix.
    // This constructor is used in combination with backSubstitute
    // to solve linear equations or invert a matrix.
    int nr = rhs.nrows();

    // Scale the matrix (find largest element of each row).
    Array1D<double> scaleVector(nr);
    for(int i = 0; i < nr; i++) {
        double maximum = 0.0;

        for(int j = 0; j < nr; j++) {
            double temp = std::abs(decomp[i][j]);
            if(temp > maximum) maximum = temp;
        }

        if(maximum == 0.0) throw SingularMatrixError("LUMatrix::LUMatrix()");
        scaleVector[i] = 1.0 / maximum;
    }

    // The loop over columns of Crout's method.
    for(int j = 0; j < nr; j++) {
        // Eqn 2.3.12 except for j = i:
        for(int i = 0; i < j; i++) {
            for(int k = 0; k < i; k++) {
                decomp[i][j] -= decomp[i][k] * decomp[k][j];
            }
        }

        // Initialize the search for the largest pivot element.
        double maximum = 0.0;
        int col_max = 0;

        // i = j of eq 2.3.12 & i = j + 1 .. N of 2.3.13.
        for(int i = j; i < nr; i++) {
            for(int k = 0; k < j; k++) {
                decomp[i][j] -= decomp[i][k] * decomp[k][j];
            }

            // Figure of merit for pivot.
            double dum = std::abs(scaleVector[i] * decomp[i][j]);

            // Is it better than the best so far ?
            if(dum >= maximum) {
                col_max = i;
                maximum = dum;
            }
        }

        // Do we need to interchange rows and the scale factor ?
        if(j != col_max)  {
            decomp.swapRows(col_max, j);
            std::swap(scaleVector[col_max], scaleVector[j]);
        }

        index[j] = col_max;

        // Now, finally, divide by the pivot element.
        if(j != nr - 1) {
            T dum = T(1) / decomp[j][j];
            for(int i = j + 1; i < nr; i++) decomp[i][j] *= dum;
        }
    }
}

template <class T>
LUMatrix<T>::~LUMatrix()
{}


template <class T>
LUMatrix<T> &LUMatrix<T>::operator=(const LUMatrix<T> &rhs) {
    decomp = rhs.decomp;
    index  = rhs.index;
    return *this;
}


template <class T> template <class I> inline
void LUMatrix<T>::backSubstitute(I iter) const {
    // Solves the set of linear equations A*X = B.
    // Here "this" is the LU-decomp of the matrix A, determined by the
    // routine LUDecompose().
    // B is input as an iterator iter pointing the the first element of the
    // vector B; the solution is returned in the same vector. This routine
    // takes into  account the possibility that B will begin with many zero
    // elements, so it is efficient for use in matrix inversion.
    // See pp 36-37 in Press & Flannery.

    int nr = decomp.nrows();
    int ii = - 1;

    for(int i = 0; i < nr; i++) {
        int ip = index[i];
        T sum = iter[ip];
        iter[ip] = iter[i];

        if(ii >= 0) {
            for(int j = ii; j < i; j++) {
                sum -= decomp[i][j] * iter[j];
            }
        } else if(sum != 0) {
            ii = i;
        }

        iter[i] = sum;
    }

    for(int i = nr; i > 0;) {
        i--;
        T sum = iter[i];

        for(int j = i + 1; j < nr; j++) {
            sum -= decomp[i][j] * iter[j];
        }

        // store a component of the solution vector X.
        iter[i] = sum / decomp[i][i];
    }
}


template <class T>
void LUMatrix<T>::backSubstitute(Vector<T> &B) const {
    if(B.size() != decomp.nrows()) {
        throw SizeError("LUMatrix::backSubstitute()", "Inconsistent dimensions.");
    }

    backSubstitute(B.begin());
}


template <class T>
void LUMatrix<T>::backSubstitute(Matrix<T> &M) const {
    if(M.nrows() != decomp.nrows()) {
        throw SizeError("LUMatrix::backSubstitute()", "Inconsistent dimensions.");
    }

    for(int column = 0; column < M.ncols(); column++) {
        backSubstitute(M.col_begin(column));
    }
}


template <class T>
Matrix<T> LUMatrix<T>::inverse() const {
    int nr = decomp.nrows();
    Matrix<T> R(nr, nr, 1.0);
    backSubstitute(R);
    return R;
}

#endif // CLASSIC_LUMatrix_HH
