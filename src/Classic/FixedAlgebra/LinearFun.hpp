#ifndef CLASSIC_LinearFun_CC
#define CLASSIC_LinearFun_CC

// ------------------------------------------------------------------------
// $RCSfile: LinearFun.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: LinearFun<T,N>
//   Representation of a linear function with values of type T,
//   and fixed number of variables N.
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/28 13:20:54 $
// $Author: mad $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/LinearFun.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "Utilities/DivideError.h"
#include "Utilities/FormatError.h"
#include <iomanip>
#include <iostream>
#include <cstring>
#ifndef IPPL_GPLUSPLUS
#include <functional>
#endif
// Implementation of template class LinearFun<T,N>.
// ------------------------------------------------------------------------

template <class T, int N>
LinearFun<T, N>::LinearFun() {
    for(int i = 0; i <= N; ++i) data[i] = 0.0;
}


template <class T, int N>
LinearFun<T, N>::LinearFun(const T &rhs) {
    data[0] = rhs;
    for(int i = 1; i <= N; ++i) data[i] = 0.0;
}


template <class T, int N>
LinearFun<T, N>::LinearFun(int rhs) {
    data[0] = T(rhs);
    for(int i = 1; i <= N; ++i) data[i] = 0.0;
}


template <class T, int N>
LinearFun<T, N> &LinearFun<T, N>::operator=(const T &rhs) {
    data[0] = rhs;
    for(int i = 1; i <= N; ++i) data[i] = 0.0;
    return *this;
}


template <class T, int N>
const T LinearFun<T, N>::getCoefficient(int index) const {
    return (index > N) ? T(0) : data[index];
}


template <class T, int N>
void LinearFun<T, N>::setCoefficient(int index, const T &value) {
    if(index <= N) data[index] = value;
}


template <class T, int N> inline
const T LinearFun<T, N>::operator[](int index) const {
    return data[index];
}


template <class T, int N> inline
T &LinearFun<T, N>::operator[](int index) {
    return data[index];
}


template <class T, int N>
LinearFun<T, N> LinearFun<T, N>::makeVariable(int var) {
    LinearFun<T, N> z;
    for(int i = 0; i <= N; ++i) z.data[i] = 0.0;
    z.data[var+1] = T(1);
    return z;
}


template <class T, int N> inline
LinearFun<T, N> LinearFun<T, N>::operator+() const {
    return *this;
}


template <class T, int N>
LinearFun<T, N> LinearFun<T, N>::operator-() const {
    LinearFun<T, N> z;
    for(int i = 0; i <= N; ++i) z.data[i] = - data[i];
    return z;
}


template <class T, int N>
LinearFun<T, N> &LinearFun<T, N>::operator+=(const LinearFun<T, N> &rhs) {
    for(int i = 0; i <= N; ++i) data[i] += rhs.data[i];
    return *this;
}


template <class T, int N>
LinearFun<T, N> &LinearFun<T, N>::operator-=(const LinearFun<T, N> &rhs) {
    for(int i = 0; i <= N; ++i) data[i] -= rhs.data[i];
    return *this;
}


template <class T, int N>
LinearFun<T, N> &LinearFun<T, N>::operator*=(const LinearFun<T, N> &rhs) {
    return *this = multiply(rhs);
}


template <class T, int N>
LinearFun<T, N> &LinearFun<T, N>::operator/=(const LinearFun<T, N> &rhs) {
    if(rhs.data[0] == T(0)) throw DivideError("LinearFun::operator/=()");
    T div = T(1) / rhs.data[0];
    data[0] *= div;

    for(int i = 1; i <= N; ++i) {
        data[i] = (data[i] - data[0] * rhs.data[i]) * div;
    }

    return *this;
}


template <class T, int N> inline
LinearFun<T, N> &LinearFun<T, N>::operator+=(const T &rhs) {
    data[0] += rhs;
    return *this;
}


template <class T, int N> inline
LinearFun<T, N> &LinearFun<T, N>::operator-=(const T &rhs) {
    data[0] -= rhs;
    return *this;
}


template <class T, int N>
LinearFun<T, N> &LinearFun<T, N>::operator*=(const T &rhs) {
    for(int i = 0; i <= N; ++i) data[i] *= rhs;
    return *this;
}


template <class T, int N>
LinearFun<T, N> &LinearFun<T, N>::operator/=(const T &rhs) {
    if(rhs == T(0)) throw DivideError("LinearFun::operator/=()");
    for(int i = 0; i <= N; ++i) data[i] /= rhs;
    return *this;
}


template <class T, int N>
bool LinearFun<T, N>::operator==(const LinearFun<T, N> &rhs) const {
    for(int i = 0; i <= N; ++i) {
        if(data[i] != rhs.data[i]) return false;
    }

    return true;
}


template <class T, int N>
bool LinearFun<T, N>::operator==(const T &rhs) const {
    if(data[0] != rhs) return false;

    for(int i = 1; i <= N; ++i) {
        if(data[i] != T(0)) return false;
    }

    return true;
}


template <class T, int N> inline
bool LinearFun<T, N>::operator!=(const LinearFun<T, N> &rhs) const {
    return !(*this == rhs);
}


template <class T, int N> inline
bool LinearFun<T, N>::operator!=(const T &rhs) const {
    return !(*this == rhs);
}


template <class T, int N>
LinearFun<T, N> LinearFun<T, N>::inverse() const {
    T aZero = data[0];
    if(aZero == T(0)) throw DivideError("LinearFun::inverse()");
    LinearFun<T, N> z;
    z.data[0] = T(1) / aZero;

    for(int i = 1; i <= N; ++i) {
        z.data[i] = - data[i] * data[i] * z.data[0];
    }

    return z;
}


template <class T, int N> LinearFun<T, N>
LinearFun<T, N>::multiply(const LinearFun<T, N> &rhs) const {
    LinearFun<T, N> z;
    z.data[0] = data[0] * rhs.data[0];

    for(int i = 1; i <= N; ++i) {
        z.data[i] = data[0] * rhs.data[i] + data[i] * rhs.data[0];
    }

    return z;
}


template <class T, int N>
T LinearFun<T, N>::evaluate(const FVector<T, N> &rhs) const {
    T z = data[0];
    for(int i = 0; i < N; ++i) z += data[i+1] * rhs[i];
    return z;
}


template <class T, int N>
LinearFun<T, N> LinearFun<T, N>::substitute(const FMatrix<T, N, N> &M) const {
    LinearFun<T, N> z;
    z.data[0] = data[0];

    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++i) {
            z.data[i+1] += data[j+1] * M[j][i];
        }
    }

    return z;
}


template <class T, int N>
LinearFun<T, N> LinearFun<T, N>::substitute(const LinearMap<T, N> &m) const {
    LinearFun<T, N> z;
    z.data[0] = data[0];

    for(int i = 0; i <= N; ++i) {
        for(int j = 0; j < N; ++i) {
            z.data[i] += data[j] * m[j][i];
        }
    }

    return z;
}


template <class T, int N>
LinearFun<T, N> LinearFun<T, N>::taylor(const T series[2]) const {
    LinearFun<T, N> z;
    z.data[0] = series[0];

    for(int i = 1; i <= N; ++i) {
        z.data[i] = series[1] * data[i];
    }

    return z;
}


template <class T, int N>
std::istream &LinearFun<T, N>::get(std::istream &is) {
    is.flags(std::ios::skipws);
    char head[4];
    is.get(head, 4);

    if(strcmp(head, "LinearFun") != 0) {
        throw FormatError("LinearFun::get()", "Flag word \"Tps\" missing.");
    }

    int maxOrder, truncOrder, nVar;
    is >> maxOrder >> truncOrder >> nVar;

    if(nVar != N) {
        throw FormatError("LinearFun::get()", "Invalid number of variables.");
    }

    LinearFun<T, N> z;
    T coeff;
    bool done = false;
    bool fail = false;

    while(true) {
        is >> coeff;
        fail = is.fail();

        int order = 0;
        int index = 0;
        for(int var = 0; var < N; var++) {
            int p;
            is >> p;
            fail |= is.fail();
            if(p < 0) {
                done = true;
            } else if(p > 0) {
                index = var;
            }

            order += p;
        }

        if(done) break;
        if(fail) throw FormatError("LinearFun::get()", "File read error");

        if(order <= 1  &&  coeff != T(0)) {
            z.data[index] = coeff;
        }
    }

    *this = z;
}


template <class T, int N>
std::ostream &LinearFun<T, N>::put(std::ostream &os) const {
    os << "LinearFun " << 1 << ' ' << 1 << ' ' << N << std::endl;
    os.precision(14);
    os.setf(std::ios::scientific, std::ios::floatfield);

    for(int i = 0; i <= N; i++) {
        if(data[i] != T(0)) {
            os << std::setw(24) << data[i];

            for(int v = 1; v <= N; ++v) {
                os << std::setw(3) << (v == i ? 1 : 0);
            }

            os << std::endl;
        }
    }

    os << std::setw(24) << T(0);

    for(int v = 0; v < N; v++) {
        os << std::setw(3) << (-1);
    }

    os << std::endl;
    os.setf(std::ios::fixed, std::ios::floatfield);
    return os;
}


// Global functions acting on LinearFun<T> objects.
// ------------------------------------------------------------------------

template <class T, int N> inline
LinearFun<T, N> operator+(const LinearFun<T, N> &lhs, const LinearFun<T, N> &rhs) {
    LinearFun<T, N> z(lhs);
    return z += rhs;
}


template <class T, int N> inline
LinearFun<T, N> operator-(const LinearFun<T, N> &lhs, const LinearFun<T, N> &rhs) {
    LinearFun<T, N> z(lhs);
    return z -= rhs;
}


template <class T, int N> inline
LinearFun<T, N> operator+(const LinearFun<T, N> &lhs, const T &rhs) {
    LinearFun<T, N> z(lhs);
    return z += rhs;
}


template <class T, int N> inline
LinearFun<T, N> operator-(const LinearFun<T, N> &lhs, const T &rhs) {
    LinearFun<T, N> z(lhs);
    return z -= rhs;
}


template <class T, int N> inline
LinearFun<T, N> operator+(const T &lhs, const LinearFun<T, N> &rhs) {
    LinearFun<T, N> z(rhs);
    return z += lhs;
}


template <class T, int N> inline
LinearFun<T, N> operator-(const T &lhs, const LinearFun<T, N> &rhs) {
    LinearFun<T, N> z(- rhs);
    return z += lhs;
}


template <class T, int N> inline
LinearFun<T, N> operator*(const LinearFun<T, N> &lhs, const LinearFun<T, N> &rhs) {
    return lhs.multiply(rhs);
}


template <class T, int N> inline
LinearFun<T, N> operator/(const LinearFun<T, N> &lhs, const LinearFun<T, N> &rhs) {
    return lhs.multiply(rhs.inverse());
}


template <class T, int N> inline
LinearFun<T, N> operator*(const LinearFun<T, N> &lhs, const T &rhs) {
    LinearFun<T, N> z(lhs);
    return z *= rhs;
}


template <class T, int N> inline
LinearFun<T, N> operator/(const LinearFun<T, N> &lhs, const T &rhs) {
    LinearFun<T, N> z(lhs);
    return z /= rhs;
}


template <class T, int N> inline
LinearFun<T, N> operator*(const T &lhs, const LinearFun<T, N> &rhs) {
    LinearFun<T, N> z(rhs);
    return z *= lhs;
}


template <class T, int N> inline
LinearFun<T, N> operator/(const T &lhs, const LinearFun<T, N> &rhs) {
    return rhs.inverse() * lhs;
}


template <class T, int N> inline
bool operator==(const T &lhs, const LinearFun<T, N> &rhs) {
    return rhs == lhs;
}


template <class T, int N> inline
bool operator!=(const T &lhs, const LinearFun<T, N> &rhs) {
    return rhs != lhs;
}


template <class T, int N> inline
std::istream &operator>>(std::istream &is, LinearFun<T, N> &tps) {
    return tps.get(is);
}


template <class T, int N> inline
std::ostream &operator<<(std::ostream &os, const LinearFun<T, N> &tps) {
    return tps.put(os);
}

#endif // CLASSIC_LinearFun_CC
