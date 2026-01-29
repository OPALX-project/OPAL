#ifndef CLASSIC_FLUMatrix_HH
#define CLASSIC_FLUMatrix_HH

// ------------------------------------------------------------------------
// $RCSfile: FLUMatrix.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FLUMatrix
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:36 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FArray1D.h"
#include "FixedAlgebra/FMatrix.h"
#include "Utilities/SingularMatrixError.h"
#include <cmath>


// Template class FLUMatrix<T,N>
// ------------------------------------------------------------------------
/// A templated representation of a LU-decomposition.
//  Constructed from a FMatrix<T,N,N>, this class remembers the triangular
//  decomposition and can return the back-substitution or the inverse.
//  When FLUMatrix is used on a complex matrix, then the statement
//  [pre]
//     include std::abs;
//  [/pre]
//  is required in the calling file, so as to inject the [tt]abs()[/tt]
//  function to the global namespace.

template <class T, int N>
class FLUMatrix {

public:

    /// Constructor.
    //  Find triangel decomposition of the matrix [b]M[/b].
    //  Throw SingularMatrixError if [b]M[/b] is singular.
    explicit FLUMatrix(const FMatrix<T, N, N> &M);

    FLUMatrix();
    FLUMatrix(const FLUMatrix<T, N> &);
    ~FLUMatrix();
    FLUMatrix<T, N> &operator=(const FLUMatrix &);

    /// Back substitution.
    //  Perform back substitution on the vector [b]B[/b].
    //  The new [b]B[/b] is the old [b]B[/b], premultiplied by the inverse.
    void backSubstitute(FVector<T, N> &B) const;

    /// Back substitution.
    //  Perform back substitution on the matrix [b]M[/b].
    //  The new [b]M[/b] is the old [b]M[/b], premultiplied by the inverse.
    template <int M>
    void backSubstitute(FMatrix<T, N, M> &MM) const;

    /// Get inverse.
    //  Construct and return the inverse.
    FMatrix<T, N, N> inverse() const;

private:

    // Back substitution on iterator.
    template <class Iterator>
    void backSubstitute(Iterator) const;

    // The LU-decomposed matrix.
    FMatrix<T, N, N> decomp;

    // The [b]index[/b] array used for LU decomposition and back substitution.
    FArray1D<int, N> index;
};


// Template function implementation.
// ------------------------------------------------------------------------

template <class T, int N>
FLUMatrix<T, N>::FLUMatrix():
    decomp(), index()
{}


template <class T, int N>
FLUMatrix<T, N>::FLUMatrix(const FLUMatrix<T, N> &rhs):
    decomp(rhs.decomp), index(rhs.index)
{}


template <class T, int N>
FLUMatrix<T, N>::FLUMatrix(const FMatrix<T, N, N> &rhs):
    decomp(rhs), index(rhs.nrows()) {
    // Builds the L-U decomp of the square matrix rhs.
    // This constructor is used in combination with backSubstitute
    // to solve linear equations or invert a matrix.

    // Scale the matrix (find largest element of each row).
    FArray1D<double, N> scaleVector;
    for(int i = 0; i < N; i++) {
        double maximum = 0.0;

        for(int j = 0; j < N; j++) {
            double temp = std::abs(decomp[i][j]);
            if(temp > maximum) maximum = temp;
        }

        if(maximum == 0.0) throw SingularMatrixError("FLUMatrix::FLUMatrix()");
        scaleVector[i] = 1.0 / maximum;
    }

    // The loop over coFLUMns of Crout's method.
    for(int j = 0; j < N; j++) {
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
        for(int i = j; i < N; i++) {
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
        if(j != N - 1) {
            T dum = T(1) / decomp[j][j];
            for(int i = j + 1; i < N; i++) decomp[i][j] *= dum;
        }
    }
}

template <class T, int N>
FLUMatrix<T, N>::~FLUMatrix()
{}


template <class T, int N>
FLUMatrix<T, N> &FLUMatrix<T, N>::operator=(const FLUMatrix<T, N> &rhs) {
    decomp = rhs.decomp;
    index  = rhs.index;
    return *this;
}


template <class T, int N> template <class I> inline
void FLUMatrix<T, N>::backSubstitute(I iter) const {
    // Solves the set of N linear equations A*X = B.
    // Here "this" is the LU-decomp of the matrix A, determined by the
    // routine LUDecompose().
    // B is input as the rhs vector B, and returns with the solution vector X.
    // This routine takes into  account the possibility that B will begin with
    // many zero elements, so it is efficient for use in matrix inversion.
    // See pp 36-37 in Press & Flannery.

    int ii = - 1;

    for(int i = 0; i < N; i++) {
        int ip = index[i];
        T sum = iter[ip];
        iter[ip] = iter[i];

        if(ii >= 0) {
            for(int j = ii; j < i; j++) {
                sum -= decomp[i][j] * iter[j];
            }
        } else if(sum != T(0)) {
            ii = i;
        }

        iter[i] = sum;
    }

    for(int i = N; i > 0;) {
        i--;
        T sum = iter[i];

        for(int j = i + 1; j < N; j++) {
            sum -= decomp[i][j] * iter[j];
        }

        // store a component of the solution vector X.
        iter[i] = sum / decomp[i][i];
    }
}


template <class T, int N>
void FLUMatrix<T, N>::backSubstitute(FVector<T, N> &B) const {
    backSubstitute(B.begin());
}


template <class T, int N> template <int M>
void FLUMatrix<T, N>::backSubstitute(FMatrix<T, N, M> &matrix) const {
    for(int column = 0; column < M; column++) {
        backSubstitute(matrix.col_begin(column));
    }
}


template <class T, int N>
FMatrix<T, N, N> FLUMatrix<T, N>::inverse() const {
    FMatrix<T, N, N> result;
    result = result + T(1);
    backSubstitute(result);
    return result;
}

#endif // CLASSIC_FLUMatrix_HH
