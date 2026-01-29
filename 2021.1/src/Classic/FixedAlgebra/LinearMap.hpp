#ifndef CLASSIC_LinearMap_CC
#define CLASSIC_LinearMap_CC

// ------------------------------------------------------------------------
// $RCSfile: LinearMap.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template Class: LinearMap<T,N>
//   Linear map in variables of type T, N variables to N variables.
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FLUMatrix.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/LinearFun.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/FVps.h"
#include "Utilities/FormatError.h"
#include "Utilities/CLRangeError.h"
#include <iostream>


// Template class LinearMap<T,N>
// ------------------------------------------------------------------------

template <class T, int N>
LinearMap<T, N>::LinearMap() {
    for(int i = 0; i < N; ++i) data[i][i+1] = 1.0;
}


template <class T, int N>
LinearMap<T, N>::LinearMap(const FVps<T, N> &rhs) {
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j <= N; ++j) data[i][j] = rhs[i][j];
    }
}


template <class T, int N>
LinearMap<T, N>::LinearMap(const FMatrix<T, N, N> &x) {
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; j++) data[i][j+1] = x[i][j];
    }
}


template <class T, int N>
LinearMap<T, N>::LinearMap(const FVector<T, N> &x) {
    for(int i = 0; i < N; i++) data[i][0] = x[i];
}


template <class T, int N>
const LinearFun<T, N> &LinearMap<T, N>::getComponent(int index) const {
    if(index > N) {
        throw CLRangeError("LinearMap::getComponent()", "Index out of range.");
    }

    return data[index];
}


template <class T, int N>
void LinearMap<T, N>::setComponent(int index, const LinearFun<T, N> &value) {
    if(index > N) {
        throw CLRangeError("LinearMap::setComponent()", "Index out of range.");
    }

    data[index] = value;
}


template <class T, int N> inline
LinearFun<T, N> &LinearMap<T, N>::operator[](int index) {
    return data[index];
}

template <class T, int N> inline
const LinearFun<T, N> &LinearMap<T, N>::operator[](int index) const {
    return data[index];
}


template <class T, int N>
LinearMap<T, N> LinearMap<T, N>::operator+() const {
    return *this;
}


template <class T, int N>
LinearMap<T, N> LinearMap<T, N>::operator-() const {
    LinearMap<T, N> z;
    for(int i = 0; i < N; i++) z[i] = - data[i];
    return z;
}


template <class T, int N>
LinearMap<T, N> &LinearMap<T, N>::operator*=(const LinearFun<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] *= rhs;
    return *this;
}


template <class T, int N>
LinearMap<T, N> &LinearMap<T, N>::operator/=(const LinearFun<T, N> &rhs) {
    LinearFun<T, N> t = rhs.inverse();
    for(int i = 0; i < N; i++) data[i] *= t;
    return *this;
}


template <class T, int N>
LinearMap<T, N> &LinearMap<T, N>::operator*=(const T &rhs) {
    for(int i = 0; i < N; i++) data[i] *= rhs;
    return *this;
}


template <class T, int N>
LinearMap<T, N> &LinearMap<T, N>::operator/=(const T &rhs) {
    for(int i = 0; i < N; i++) data[i] /= rhs;
    return *this;
}


template <class T, int N>
LinearMap<T, N>& LinearMap<T, N>::operator+=(const LinearMap<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] += rhs[i];
    return *this;
}


template <class T, int N>
LinearMap<T, N>& LinearMap<T, N>::operator-=(const LinearMap<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] -= rhs[i];
    return *this;
}


template <class T, int N>
LinearMap<T, N>& LinearMap<T, N>::operator+=(const FVector<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] += rhs[i];
    return *this;
}


template <class T, int N>
LinearMap<T, N>& LinearMap<T, N>::operator-=(const FVector<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] -= rhs[i];
    return *this;
}


template <class T, int N>
std::istream &LinearMap<T, N>::get(std::istream &is) {
    is.flags(std::ios::skipws);
    char head[4];
    (is >> std::ws).get(head, 4);

    if(strcmp(head, "LinearMap") != 0) {
        throw FormatError("LinearMap::get()", "Flag word \"LinearMap\" missing.");
    }

    int nDim;
    is >> nDim;

    if(nDim != N) {
        throw FormatError("LinearMap::get()", "Invalid LinearMap dimension");
    }

    // Read into temporary for exception safety.
    LinearMap z;
    for(int i = 0; i < N; i++) is >> z.data[i];
    *this = z;
    return is;
}


template <class T, int N>
std::ostream &LinearMap<T, N>::put(std::ostream &os) const {
    os << "LinearMap " << N << std::endl;
    for(int i = 0; i < N; i++) os << data[i];
    return os;
}


template <class T, int N>
LinearMap<T, N> LinearMap<T, N>::inverse() const {
    FVector<T, N> vec = constantTerm();
    FLUMatrix<T, N> lu(linearTerms());
    FMatrix<T, N, N> mat(lu.inverse());

    // Initialize map z with inverse linear map.
    LinearMap<T, N> z(mat);
    z -= mat * vec;
    return z;
}


template <class T, int N>
void LinearMap<T, N>::identity() {
    for(int i = 0; i < N; ++i) data[i] = LinearFun<T, N>::makeVariable(i);
}


template <class T, int N>
FVector<T, N> LinearMap<T, N>::constantTerm(const FVector<T, N> &rhs) const {
    FVector<T, N> z;

    for(int i = 0; i < N; ++i) {
        z[i] = data[i][0];
        for(int j = 0; j < N; ++j) {
            z[i] += data[i][j+1] * rhs[j];
        }
    }

    return z;
}


template <class T, int N>
FVector<T, N> LinearMap<T, N>::constantTerm() const {
    FVector<T, N> z;
    for(int i = 0; i < N; i++) z[i] = data[i][0];
    return z;
}


template <class T, int N>
FMatrix<T, N, N> LinearMap<T, N>::linearTerms() const {
    FMatrix<T, N, N> z;

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            z(i, j) = data[i][j+1];
        }
    }

    return z;
}


template <class T, int N>
LinearMap<T, N> LinearMap<T, N>::substitute(const LinearMap<T, N> &rhs) const {
    LinearMap<T, N> z;

    for(int i = 0; i < N; ++i) {
        z.data[i][0] = data[i][0];
        for(int j = 0; j < N; ++j) {
            for(int k = 0; k <= N; ++k) {
                z.data[i][k] += data[i][j+1] * rhs.data[j][k];
            }
        }
    }

    return z;
}


template <class T, int N>
LinearMap<T, N> LinearMap<T, N>::substitute(const FMatrix<T, N, N> &rhs) const {
    LinearMap<T, N> z;

    for(int i = 0; i < N; ++i) {
        z.data[i][0] = data[i][0];
        for(int j = 0; j < N; ++j) {
            for(int k = 0; k < N; ++k) {
                z.data[i][k+1] += data[i][j+1] * rhs[j][k];
            }
        }
    }

    return z;
}


template <class T, int N>
LinearMap<T, N> LinearMap<T, N>::substituteInto(const FMatrix<T, N, N> &lhs) const {
    LinearMap<T, N> z;

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            z[i] += lhs[i][j] * data[j];
        }
    }

    return z;
}


// Global Operators on LinearMap<T,N>
// ------------------------------------------------------------------------

template <class T, int N>
LinearMap<T, N> operator*(const LinearMap<T, N> &lhs, const LinearFun<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] * rhs;
    return z;
}


template <class T, int N>
LinearMap<T, N> operator*(const LinearFun<T, N> &lhs, const LinearMap<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs * rhs[i];
    return z;
}


template <class T, int N>
LinearMap<T, N> operator*(const LinearMap<T, N> &lhs, const T &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] * rhs;
    return z;
}


template <class T, int N>
LinearMap<T, N> operator*(const T &lhs, const LinearMap<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs * rhs[i];
    return z;
}

template <class T, int N>
LinearMap<T, N> operator*(const FMatrix<T, N, N> &lhs, const LinearMap<T, N> &rhs) {
    return rhs.substituteInto(lhs);
}


template <class T, int N>
LinearMap<T, N> operator/(const LinearMap<T, N> &lhs, const LinearFun<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] / rhs;
    return z;
}


template <class T, int N>
LinearMap<T, N> operator/(const LinearMap<T, N> &lhs, const T &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] / rhs;
    return z;
}


template <class T, int N>
LinearMap<T, N> operator+(const LinearMap<T, N> &lhs, const LinearMap<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] + rhs[i];
    return z;
}


template <class T, int N>
LinearMap<T, N> operator-(const LinearMap<T, N> &lhs, const LinearMap<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] - rhs[i];
    return z;
}


template <class T, int N>
LinearMap<T, N> operator+(const LinearMap<T, N> &lhs, const FVector<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] + rhs[i];
    return z;
}


template <class T, int N>
LinearMap<T, N> operator-(const LinearMap<T, N> &lhs, const FVector<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] - rhs[i];
    return z;
}


template <class T, int N>
LinearMap<T, N> operator+(const FVector<T, N> &lhs, const LinearMap<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] + rhs[i];
    return z;
}


template <class T, int N>
LinearMap<T, N> operator-(const FVector<T, N> &lhs, const LinearMap<T, N> &rhs) {
    LinearMap<T, N> z;
    for(int i = 0; i < N; ++i) z[i] = lhs[i] - rhs[i];
    return z;
}


template <class T, int N>
std::istream &operator>>(std::istream &is, LinearMap<T, N> &vps) {
    return vps.get(is);
}


template <class T, int N>
std::ostream &operator<<(std::ostream &os, const LinearMap<T, N> &vps) {
    return vps.put(os);
}

#endif // CLASSIC_LinearMap_CC
