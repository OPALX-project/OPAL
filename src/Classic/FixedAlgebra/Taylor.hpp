#ifndef CLASSIC_Taylor_CC
#define CLASSIC_Taylor_CC
// ------------------------------------------------------------------------
// $RCSfile: Taylor.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: Taylor<T>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include <numeric>
#include <iomanip>

#include "FixedAlgebra/Taylor.h"


// Template class Taylor<T>.
// ------------------------------------------------------------------------

template <class T>
Taylor<T>::
Taylor(int order):
    itsCoeffs(order + 1)
{}


template <class T>
Taylor<T>::
Taylor():
    itsCoeffs()
{}


template <class T>
Taylor<T>::
Taylor(const Taylor &rhs):
    itsCoeffs(rhs.itsCoeffs)
{}


template <class T>
Taylor<T>::
~Taylor()
{}


template <class T>
const Taylor<T> &Taylor<T>::
operator=(const Taylor &rhs) {
    itsCoeffs = rhs.itsCoeffs;
    return *this;
}


template<class T> inline
T *Taylor<T>::
begin() {
    return itsCoeffs.begin();
}


template<class T> inline
const T *Taylor<T>::
begin() const {
    return itsCoeffs.begin();
}


template<class T> inline
T *Taylor<T>::
end() {
    return itsCoeffs.end();
}


template<class T> inline
const T *Taylor<T>::
end() const {
    return itsCoeffs.end();
}


template<class T> inline
T &Taylor<T>::operator[](int i) {
    return itsCoeffs[i];
}


template<class T> inline
const T &Taylor<T>::operator[](int i) const {
    return itsCoeffs[i];
}


template<class T>
Taylor<T> Taylor<T>::operator-() const {
    Taylor<T> result(getOrder());
    std::transform(begin(), end(), result.begin(), std::negate<T>());
    return result;
}


template<class T>
Taylor<T> &Taylor<T>::operator*=(const T &val) {
    std::transform(begin(), end(), begin(),
                   std::bind2nd(std::multiplies<T>(), val));
    return *this;
}


template<class T>
Taylor<T> &Taylor<T>::operator/=(const T &val) {
    std::transform(begin(), end(), begin(),
                   std::bind2nd(std::divides<T>(), val));
    return *this;
}


template<class T>
Taylor<T> &Taylor<T>::operator+=(const Taylor<T> &rhs) {
    if(getOrder() < rhs.getOrder()) {
        Taylor<T> temp(rhs.getOrder());
        std::copy(begin(), end(), temp.begin());
        *this = temp;
    }

    int n = rhs.end() - rhs.begin();
    std::transform(begin(), begin() + n, rhs.begin(), begin(), std::plus<T>());
    return *this;
}


template<class T>
Taylor<T> &Taylor<T>::operator-=(const Taylor<T> &rhs) {
    if(getOrder() < rhs.getOrder()) {
        Taylor<T> temp(rhs.getOrder());
        std::copy(begin(), end(), temp.begin());
        *this = temp;
    }

    int n = rhs.end() - rhs.begin();
    std::transform(begin(), begin() + n, rhs.begin(), begin(), std::minus<T>());
    return *this;
}


template <class T>
void Taylor<T>::
clear() {
    erase(begin(), end());
}


template <class T>
Taylor<T> Taylor<T>::
integrate() const {
    Taylor<T> result(getOrder() + 1);
    const T *x1 = begin();
    const T *x2 = end();
    T *z = result.begin();
    double div = 0.0;

    for(const T *x = x1; x < x2; ++x, ++z) {
        div += 1.0;
        *z = *x / div;
    }

    return result;
}


template <class T>
T Taylor<T>::
sum() const {
    return std::accumulate(begin(), end(), T(0));
}


template <class T> inline
int Taylor<T>::
getOrder() const {
    return itsCoeffs.size() - 1;
}


// Global functions.
// ------------------------------------------------------------------------

template <class T>
Taylor<T>
operator+(const Taylor<T> &lhs, const Taylor<T> &rhs) {
    Taylor<T> result(std::max(lhs.getOrder(), rhs.getOrder()));
    const T *x = lhs.begin();
    const T *y = rhs.begin();
    T *z = result.begin();
    while(x != lhs.end()  &&  y != rhs.end()) {
        *z++ = *x++ + *y++;
    }
    while(x != lhs.end()) *z++ = *x++;
    while(y != rhs.end()) *z++ = *y++;
    return result;
}


template <class T>
Taylor<T>
operator-(const Taylor<T> &lhs, const Taylor<T> &rhs) {
    Taylor<T> result(std::max(lhs.getOrder(), rhs.getOrder()));
    const T *x = lhs.begin();
    const T *y = rhs.begin();
    T *z = result.begin();
    while(x != lhs.end()  &&  y != rhs.end()) {
        *z++ = *x++ - *y++;
    }
    while(x != lhs.end()) *z++ = *x++;
    while(y != rhs.end()) *z++ = - *y++;
    return result;
}


template <class T>
Taylor<T>
operator*(const Taylor<T> &lhs, double rhs) {
    Taylor<T> result(lhs.getOrder());
    const T *x = lhs.begin();
    T *z = result.begin();
    while(x != lhs.end()) {
        *z++ = *x++ * rhs;
    }
    return result;
}


template <class T>
Taylor<T>
operator*(double lhs, const Taylor<T> &rhs) {
    Taylor<T> result(rhs.getOrder());
    const T *x = rhs.begin();
    T *z = result.begin();
    while(x != rhs.end()) {
        *z++ = lhs * *x++;
    }
    return result;
}


template <class T>
Taylor<T>
operator/(const Taylor<T> &lhs, double rhs) {
    Taylor<T> result(lhs.getOrder());
    const T *x = lhs.begin();
    T *z = result.begin();
    while(x != lhs.end()) {
        *z++ = *x++ / rhs;
    }
    return result;
}


template <class T>
Taylor<T>
PoissonBracket(const Taylor<T> &lhs, const Taylor<T> &rhs) {
    Taylor<T> result(lhs.getOrder() + rhs.getOrder());
    for(int i = 0; i <= lhs.getOrder(); ++i) {
        if(! lhs[i].isZero()) {
            for(int j = 0; j <= rhs.getOrder(); ++j) {
                if(! rhs[j].isZero()) {
                    result[i+j] += PoissonBracket(lhs[i], rhs[j]);
                }
            }
        }
    }
    return result;
}


template <class T>
std::ostream &operator<<(std::ostream &os, const Taylor<T> &series) {
    os << "Taylor" << std::setw(4) << series.getOrder() << std::endl;

    for(int i = 0; i <= series.getOrder(); ++i) {
        os << "Terms of order" << std::setw(4) << i << " " << series[i];
    }

    return os;
}

#endif
