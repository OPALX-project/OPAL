#ifndef CLASSIC_TransportFun_CC
#define CLASSIC_TransportFun_CC

// ------------------------------------------------------------------------
// $RCSfile: TransportFun.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: TransportFun<T,N>
//   Representation of a second-order function with values of type T,
//   and fixed number of variables N.
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/TransportFun.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/TransportMap.h"
#include "Utilities/DivideError.h"
#include "Utilities/FormatError.h"
#include <iomanip>
#include <iostream>
#include <cstring>


// Implementation of template class TransportFun<T,N>.
// ------------------------------------------------------------------------

template <class T, int N>
TransportFun<T, N>::TransportFun() {
    for(int i = 0; i < SIZE; ++i) data[i] = 0.0;
}


template <class T, int N>
TransportFun<T, N>::TransportFun(const T &rhs) {
    data[0] = rhs;
    for(int i = 1; i < SIZE; ++i) data[i] = 0.0;
}


template <class T, int N>
TransportFun<T, N>::TransportFun(int rhs) {
    data[0] = T(rhs);
    for(int i = 1; i < SIZE; ++i) data[i] = 0.0;
}


template <class T, int N>
TransportFun<T, N> &TransportFun<T, N>::operator=(const T &rhs) {
    data[0] = rhs;
    for(int i = 1; i < SIZE; ++i) data[i] = 0.0;
    return *this;
}


template <class T, int N>
const T TransportFun<T, N>::getCoefficient(int index) const {
    return (index > N) ? T(0) : data[index];
}


template <class T, int N>
const T TransportFun<T, N>::getCoefficient(int index1, int index2) const {
    if(index1 <= index2  &&  index2 < N) {
        int index = N * (2 * N + 1 - index1) * index1 + (N + 1) + index2;
        return data[index];
    } else {
        return T(0);
    }
}


template <class T, int N>
void TransportFun<T, N>::setCoefficient(int index, const T &value) {
    if(index <= N) data[index] = value;
}


template <class T, int N>
void TransportFun<T, N>::setCoefficient(int index1, int index2, const T &value) {
    if(index1 <= index2  &&  index2 < N) {
        int index = N * (2 * N + 1 - index1) * index1 + (N + 1) + index2;
        data[index] = value;
    }
}


template <class T, int N> inline
const T TransportFun<T, N>::operator[](int index) const {
    return data[index];
}


template <class T, int N> inline
T &TransportFun<T, N>::operator[](int index) {
    return data[index];
}


template <class T, int N> inline
const T TransportFun<T, N>::operator()(int index1, int index2) const {
    int index = N * (2 * N + 1 - index1) * index1 + (N + 1) + index2;
    return data[index];
}


template <class T, int N> inline
T &TransportFun<T, N>::operator()(int index1, int index2) {
    int index = N * (2 * N + 1 - index1) * index1 + (N + 1) + index2;
    return data[index];
}


template <class T, int N>
TransportFun<T, N> TransportFun<T, N>::makeVariable(int var) {
    TransportFun<T, N> z;
    for(int i = 0; i < SIZE; ++i) z.data[i] = 0.0;
    z.data[var+1] = T(1);
    return z;
}


template <class T, int N> inline
TransportFun<T, N> TransportFun<T, N>::operator+() const {
    return *this;
}


template <class T, int N>
TransportFun<T, N> TransportFun<T, N>::operator-() const {
    TransportFun<T, N> z;
    for(int i = 0; i < SIZE; ++i) z.data[i] = - data[i];
    return z;
}


template <class T, int N>
TransportFun<T, N> &TransportFun<T, N>::operator+=(const TransportFun<T, N> &rhs) {
    for(int i = 0; i < SIZE; ++i) data[i] += rhs.data[i];
    return *this;
}


template <class T, int N>
TransportFun<T, N> &TransportFun<T, N>::operator-=(const TransportFun<T, N> &rhs) {
    for(int i = 0; i < SIZE; ++i) data[i] -= rhs.data[i];
    return *this;
}


template <class T, int N>
TransportFun<T, N> &TransportFun<T, N>::operator*=(const TransportFun<T, N> &rhs) {
    return *this = multiply(rhs);
}


template <class T, int N>
TransportFun<T, N> &TransportFun<T, N>::operator/=(const TransportFun<T, N> &rhs) {
    if(rhs.data[0] == T(0)) throw DivideError("TransportFun::operator/=()");
    return *this = multiply(rhs.inverse());
}


template <class T, int N> inline
TransportFun<T, N> &TransportFun<T, N>::operator+=(const T &rhs) {
    data[0] += rhs;
    return *this;
}


template <class T, int N> inline
TransportFun<T, N> &TransportFun<T, N>::operator-=(const T &rhs) {
    data[0] -= rhs;
    return *this;
}


template <class T, int N>
TransportFun<T, N> &TransportFun<T, N>::operator*=(const T &rhs) {
    for(int i = 0; i < SIZE; ++i) data[i] *= rhs;
    return *this;
}


template <class T, int N>
TransportFun<T, N> &TransportFun<T, N>::operator/=(const T &rhs) {
    if(rhs == T(0)) throw DivideError("TransportFun::operator/=()");
    for(int i = 0; i < SIZE; ++i) data[i] /= rhs;
    return *this;
}


template <class T, int N>
bool TransportFun<T, N>::operator==(const TransportFun<T, N> &rhs) const {
    for(int i = 0; i < SIZE; ++i) {
        if(data[i] != rhs.data[i]) return false;
    }

    return true;
}


template <class T, int N>
bool TransportFun<T, N>::operator==(const T &rhs) const {
    if(data[0] != rhs) return false;

    for(int i = 1; i < SIZE; ++i) {
        if(data[i] != T(0)) return false;
    }

    return true;
}


template <class T, int N> inline
bool TransportFun<T, N>::operator!=(const TransportFun<T, N> &rhs) const {
    return !(*this == rhs);
}


template <class T, int N> inline
bool TransportFun<T, N>::operator!=(const T &rhs) const {
    return !(*this == rhs);
}


template <class T, int N>
TransportFun<T, N> TransportFun<T, N>::inverse() const {
    T aZero = data[0];
    if(aZero == T(0)) throw DivideError("TransportFun::inverse()");

    T series[3];
    series[0] = T(1) / aZero;
    series[1] = - series[0] / aZero;
    series[2] = - series[1] / aZero;
    return taylor(series);
}


template <class T, int N> TransportFun<T, N>
TransportFun<T, N>::multiply(const TransportFun<T, N> &rhs) const {
    TransportFun<T, N> z;

    // Product of constant terms.
    z.data[0] = data[0] * rhs.data[0];

    // Products of constant term with non-constant terms of other factor.
    for(int i = 1; i < SIZE; ++i) {
        z.data[i] = data[0] * rhs.data[i] + data[i] * rhs.data[0];
    }

    // Second-order terms generated as products of linear terms.
    int k = N;
    for(int i = 0; i < N; ++i) {
        // Quadratic terms.
        ++k;
        z.data[k] += data[i+1] * rhs.data[i+1];

        // Bilinear terms.
        for(int j = i + 1; j < N; ++j) {
            ++k;
            z.data[k] += data[i+1] * rhs.data[j+1] + data[j+1] * rhs.data[i+1];
        }
    }

    return z;
}


template <class T, int N>
T TransportFun<T, N>::evaluate(const FVector<T, N> &rhs) const {
    // Constant term.
    T z = data[0];

    // Linear terms.
    int k = 0;
    for(int i = 0; i < N; ++i) {
        z += data[++k] * rhs[i];
    }

    // Second-order terms.
    for(int i = 0; i < N; ++i) {
        for(int j = i; j < N; ++j) {
            z += data[++k] * rhs[i] * rhs[j];
        }
    }

    return z;
}


template <class T, int N>
TransportFun<T, N> TransportFun<T, N>::substitute(const FMatrix<T, N, N> &M) const {
    return substitute(TransportMap<T, N>(M));
}


template <class T, int N>
TransportFun<T, N> TransportFun<T, N>::substitute(const TransportMap<T, N> &m) const {
    TransportFun<T, N> z;
    z.data[0] = data[0];

    // Linear terms.
    int k = 0;
    for(int i = 0; i < N; ++i) {
        z += data[++k] * m[i];
    }

    // Second-order terms.
    for(int i = 0; i < N; ++i) {
        for(int j = i; j < N; ++i) {
            z += data[++k] * m[i] * m[j];
        }
    }

    return z;
}


template <class T, int N>
TransportFun<T, N> TransportFun<T, N>::taylor(const T series[3]) const {
    TransportFun<T, N> x(*this);
    x[0] = 0.0;
    return (series[2] * x + series[1]) * x + series[0];
}


template <class T, int N>
std::istream &TransportFun<T, N>::get(std::istream &is) {
    is.flags(std::ios::skipws);
    char head[4];
    is.get(head, 4);

    if(strcmp(head, "Tps") != 0) {
        throw FormatError("TransportFun::get()", "Flag word \"Tps\" missing.");
    }

    int maxOrder, truncOrder, nVar;
    is >> maxOrder >> truncOrder >> nVar;

    if(nVar != N) {
        throw FormatError("TransportFun::get()", "Invalid number of variables.");
    }

    T coeff;
    int mono[N];
    bool done = false;
    bool fail = false;

    while(true) {
        is >> coeff;
        fail = is.fail();

        int order = 0;
        for(int var = 0; var < nVar; ++var) {
            is >> mono[var];
            fail |= is.fail();
            if(mono[var] < 0) done = true;
            order += mono[var];
        }

        if(done) break;
        if(fail) throw FormatError("TransportFun::get()", "File read error");

        // Store coefficient in proper place.
        if(order <= 2  &&  coeff != T(0)) {
            int order = 0;
            int index = 0;
            for(int var = N; var-- > 0;) {
                order += mono[var];
                if(order == 1) {
                    ++index;
                } else if(order == 2) {
                    index += N + order - var;
                }
            }
            data[index] = coeff;
        }
    }
}


template <class T, int N>
std::ostream &TransportFun<T, N>::put(std::ostream &os) const {
    os << "Tps " << this->itsRep->max << ' ' << this->itsRep->trc << ' ' << N
       << std::endl;
    std::streamsize old_prec = os.precision(14);
    os.setf(std::ios::scientific, std::ios::floatfield);

    // Constant term.
    T *coeff = data;
    if(*coeff != T(0)) {
        os << std::setw(24) << *coeff[0];
        for(int var = 0; var < N; ++var) {
            os << std::setw(3) << 0;
        }
    }

    // Linear terms.
    for(int i = 0; i < N; i++) {
        ++coeff;
        if(*coeff != T(0)) {
            os << std::setw(24) << *coeff;
            for(int var = 0; var < N; ++var) {
                os << std::setw(3) << ((var == i) ? 1 : 0);
            }
            os << std::endl;
        }
    }

    // Second-order terms.
    for(int i = 0; i < N; ++i) {
        ++coeff;
        if(*coeff != T(0)) {
            os << std::setw(24) << *coeff;
            for(int var = 0; var < N; ++var) {
                os << std::setw(3) << ((var == i) ? 2 : 0);
            }
        }
        for(int j = i + 1; j < N; ++j) {
            ++coeff;
            if(*coeff != T(0)) {
                os << std::setw(24) << *coeff;
                for(int var = 0; var < N; ++var) {
                    os << std::setw(3) << ((var == i || var == j) ? 1 : 0);
                }
                os << std::endl;
            }
        }
    }

    // End marker.
    os << std::setw(24) << T(0);
    for(int var = 0; var < N; ++var) {
        os << std::setw(3) << (-1);
    }
    os << std::endl;

    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
    return os;
}


// Global functions acting on TransportFun<T> objects.
// ------------------------------------------------------------------------

template <class T, int N> inline
TransportFun<T, N> operator+(const TransportFun<T, N> &lhs, const TransportFun<T, N> &rhs) {
    TransportFun<T, N> z(lhs);
    return z += rhs;
}


template <class T, int N> inline
TransportFun<T, N> operator-(const TransportFun<T, N> &lhs, const TransportFun<T, N> &rhs) {
    TransportFun<T, N> z(lhs);
    return z -= rhs;
}


template <class T, int N> inline
TransportFun<T, N> operator+(const TransportFun<T, N> &lhs, const T &rhs) {
    TransportFun<T, N> z(lhs);
    return z += rhs;
}


template <class T, int N> inline
TransportFun<T, N> operator-(const TransportFun<T, N> &lhs, const T &rhs) {
    TransportFun<T, N> z(lhs);
    return z -= rhs;
}


template <class T, int N> inline
TransportFun<T, N> operator+(const T &lhs, const TransportFun<T, N> &rhs) {
    TransportFun<T, N> z(rhs);
    return z += lhs;
}


template <class T, int N> inline
TransportFun<T, N> operator-(const T &lhs, const TransportFun<T, N> &rhs) {
    TransportFun<T, N> z(- rhs);
    return z += lhs;
}


template <class T, int N> inline
TransportFun<T, N> operator*(const TransportFun<T, N> &lhs, const TransportFun<T, N> &rhs) {
    return lhs.multiply(rhs);
}


template <class T, int N> inline
TransportFun<T, N> operator/(const TransportFun<T, N> &lhs, const TransportFun<T, N> &rhs) {
    return lhs.multiply(rhs.inverse());
}


template <class T, int N> inline
TransportFun<T, N> operator*(const TransportFun<T, N> &lhs, const T &rhs) {
    TransportFun<T, N> z(lhs);
    return z *= rhs;
}


template <class T, int N> inline
TransportFun<T, N> operator/(const TransportFun<T, N> &lhs, const T &rhs) {
    TransportFun<T, N> z(lhs);
    return z /= rhs;
}


template <class T, int N> inline
TransportFun<T, N> operator*(const T &lhs, const TransportFun<T, N> &rhs) {
    TransportFun<T, N> z(rhs);
    return z *= lhs;
}


template <class T, int N> inline
TransportFun<T, N> operator/(const T &lhs, const TransportFun<T, N> &rhs) {
    return rhs.inverse() * lhs;
}


template <class T, int N> inline
bool operator==(const T &lhs, const TransportFun<T, N> &rhs) {
    return rhs == lhs;
}


template <class T, int N> inline
bool operator!=(const T &lhs, const TransportFun<T, N> &rhs) {
    return rhs != lhs;
}


template <class T, int N>
std::istream &operator>>(std::istream &is, TransportFun<T, N> &fun) {
    return fun.get(is);
}


template <class T, int N>
std::ostream &operator<<(std::ostream &os, const TransportFun<T, N> &fun) {
    return fun.put(os);
}

#endif // CLASSIC_TransportFun_CC
