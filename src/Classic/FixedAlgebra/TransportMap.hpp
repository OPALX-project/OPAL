#ifndef CLASSIC_TransportMap_CC
#define CLASSIC_TransportMap_CC

// ------------------------------------------------------------------------
// $RCSfile: TransportMap.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template Class: TransportMap<T,N>
//   Second-order map in variables of type T, N variables to N variables.
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/28 13:20:56 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FLUMatrix.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/TransportFun.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/FVps.h"
#include "Utilities/FormatError.h"
#include "Utilities/CLRangeError.h"
#include <iostream>


// Template class TransportMap<T,N>
// ------------------------------------------------------------------------

template <class T, int N>
const int TransportMap<T, N>::SIZE;


template <class T, int N>
TransportMap<T, N>::TransportMap() {
    for(int i = 0; i < N; ++i) data[i][i+1] = 1.0;
}


template <class T, int N>
TransportMap<T, N>::TransportMap(const FVps<T, N> &rhs) {
    for(int i = 0; i < N; ++i) {
        int size = std::min(rhs[i].getSize(), SIZE);
        for(int j = 0; j < size; ++j) data[i][j] = rhs[i][j];
    }
}


template <class T, int N>
TransportMap<T, N>::TransportMap(const FMatrix<T, N, N> &x) {
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; j++) data[i][j+1] = x[i][j];
    }
}


template <class T, int N>
TransportMap<T, N>::TransportMap(const FVector<T, N> &x) {
    for(int i = 0; i < N; i++) data[i][0] = x[i];
}


template <class T, int N>
const TransportFun<T, N> &TransportMap<T, N>::getComponent(int index) const {
    if(index > N) {
        throw CLRangeError("TransportMap::getComponent()", "Index out of range.");
    }

    return data[index];
}


template <class T, int N>
void TransportMap<T, N>::setComponent(int index, const TransportFun<T, N> &value) {
    if(index > N) {
        throw CLRangeError("TransportMap::setComponent()", "Index out of range.");
    }

    data[index] = value;
}


template <class T, int N> inline
TransportFun<T, N> &TransportMap<T, N>::operator[](int index) {
    return data[index];
}

template <class T, int N> inline
const TransportFun<T, N> &TransportMap<T, N>::operator[](int index) const {
    return data[index];
}


template <class T, int N>
TransportMap<T, N> TransportMap<T, N>::operator+() const {
    return *this;
}


template <class T, int N>
TransportMap<T, N> TransportMap<T, N>::operator-() const {
    TransportMap<T, N> z;
    for(int i = 0; i < N; i++) z[i] = - data[i];
    return z;
}


template <class T, int N>
TransportMap<T, N> &TransportMap<T, N>::operator*=(const TransportFun<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] *= rhs;
    return *this;
}


template <class T, int N>
TransportMap<T, N> &TransportMap<T, N>::operator/=(const TransportFun<T, N> &rhs) {
    TransportFun<T, N> t = rhs.inverse();
    for(int i = 0; i < N; i++) data[i] *= t;
    return *this;
}


template <class T, int N>
TransportMap<T, N> &TransportMap<T, N>::operator*=(const T &rhs) {
    for(int i = 0; i < N; i++) data[i] *= rhs;
    return *this;
}


template <class T, int N>
TransportMap<T, N> &TransportMap<T, N>::operator/=(const T &rhs) {
    for(int i = 0; i < N; i++) data[i] /= rhs;
    return *this;
}


template <class T, int N>
TransportMap<T, N>& TransportMap<T, N>::operator+=(const TransportMap<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] += rhs[i];
    return *this;
}


template <class T, int N>
TransportMap<T, N>& TransportMap<T, N>::operator-=(const TransportMap<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] -= rhs[i];
    return *this;
}


template <class T, int N>
TransportMap<T, N>& TransportMap<T, N>::operator+=(const FVector<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] += rhs[i];
    return *this;
}


template <class T, int N>
TransportMap<T, N>& TransportMap<T, N>::operator-=(const FVector<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] -= rhs[i];
    return *this;
}


template <class T, int N>
TransportMap<T, N> TransportMap<T, N>::inverse() const {
    // Separate out const, linear, and second-order terms.
    //    C := A(0), the constant part of map A,
    //    M := A(1), the linear part of map A,
    //    S := - A(2) = - (A - M - C).
    FVector<T, N> C;
    FMatrix<T, N, N> M;
    TransportMap<T, N> S;

    for(int i = 0; i < N; ++i) {
        const TransportFun<T, N> &fun = data[i];
        C[i] = - fun[0];
        for(int j = 0; j < N; ++j) {
            M[i][j] = fun[j+1];
        }
        for(int j = N + 1; j < SIZE; ++j) {
            S[i][j] = - fun[j];
        }
    }

    FLUMatrix<T, N> lu(M);
    FMatrix<T, N, N> mat(lu.inverse());

    // Initialize map result with inverse linear map.
    TransportMap<T, N> result(C);

    // Compute second-order part of inverse.
    TransportMap<T, N> L;
    L -= C;
    result = mat * (S.substitute(mat * L) + L);

    return result;
}


template <class T, int N>
void TransportMap<T, N>::identity() {
    for(int i = 0; i < N; ++i) data[i] = TransportFun<T, N>::makeVariable(i);
}


template <class T, int N>
FVector<T, N> TransportMap<T, N>::constantTerm(const FVector<T, N> &rhs) const {
    FVector<T, N> z;

    for(int i = 0; i < N; ++i) {
        z[i] = data[i].evaluate(rhs);
    }

    return z;
}


template <class T, int N>
FVector<T, N> TransportMap<T, N>::constantTerm() const {
    FVector<T, N> z;

    for(int i = 0; i < N; i++) {
        z[i] = data[i][0];
    }

    return z;
}


template <class T, int N>
FMatrix<T, N, N> TransportMap<T, N>::linearTerms() const {
    FMatrix<T, N, N> z;

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            z(i, j) = data[i][j+1];
        }
    }

    return z;
}


template <class T, int N>
TransportMap<T, N> TransportMap<T, N>::
substitute(const TransportMap<T, N> &rhs) const {
    TransportMap<T, N> z;

    for(int n = 0; n < N; ++n) {
        // Constant terms.
        z.data[n] = data[n][0];
    }

    // Linear terms.
    for(int i = 0; i < N; ++i) {
        const TransportFun<T, N> &x = rhs[i];
        for(int n = 0; n < N; ++n) {
            z[n] += data[n][i+1] * x;
        }
    }

    // Second-order terms.
    int k = N;
    for(int i = 0; i < N; ++i) {
        for(int j = i; j < N; ++j) {
            const TransportFun<T, N> x = rhs[i] * rhs[j];
            ++k;
            for(int n = 0; n < N; ++n) {
                z[n] += data[n][k] * x;
            }
        }
    }

    return z;
}


template <class T, int N>
TransportMap<T, N> TransportMap<T, N>::
substitute(const FMatrix<T, N, N> &rhs) const {
    return substitute(TransportMap<T, N>(rhs));
}


template <class T, int N>
TransportMap<T, N> TransportMap<T, N>::
substituteInto(const FMatrix<T, N, N> &lhs) const {
    TransportMap<T, N> z;

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            z[i] += lhs[i][j] * data[j];
        }
    }

    return z;
}


// Global Operators on TransportMap<T,N>
// ------------------------------------------------------------------------

template <class T, int N>
TransportMap<T, N>
operator*(const TransportMap<T, N> &lhs, const TransportFun<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] * rhs;
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator*(const TransportFun<T, N> &lhs, const TransportMap<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs * rhs[i];
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator*(const TransportMap<T, N> &lhs, const T &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] * rhs;
    return z;
}


template <class T, int N>
TransportMap<T, N> operator*(const T &lhs, const TransportMap<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs * rhs[i];
    return z;
}

template <class T, int N>
TransportMap<T, N> operator*(const FMatrix<T, N, N> &lhs, const TransportMap<T, N> &rhs) {
    return rhs.substituteInto(lhs);
}


template <class T, int N>
TransportMap<T, N>
operator/(const TransportMap<T, N> &lhs, const TransportFun<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] / rhs;
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator/(const TransportMap<T, N> &lhs, const T &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] / rhs;
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator+(const TransportMap<T, N> &lhs, const TransportMap<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] + rhs[i];
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator-(const TransportMap<T, N> &lhs, const TransportMap<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] - rhs[i];
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator+(const TransportMap<T, N> &lhs, const FVector<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] + rhs[i];
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator-(const TransportMap<T, N> &lhs, const FVector<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] - rhs[i];
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator+(const FVector<T, N> &lhs, const TransportMap<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] + rhs[i];
    return z;
}


template <class T, int N>
TransportMap<T, N>
operator-(const FVector<T, N> &lhs, const TransportMap<T, N> &rhs) {
    TransportMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] - rhs[i];
    return z;
}

#endif // CLASSIC_TransportMap_CC
