#ifndef CLASSIC_FLieGenerator_CC
#define CLASSIC_FLieGenerator_CC
// ------------------------------------------------------------------------
// $RCSfile: FLieGenerator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FLieGenerator<T,N>
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2001/11/15 08:52:26 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "FixedAlgebra/FLieGenerator.h"
#include "Algebra/Array1D.h"
#include "FixedAlgebra/FArray1D.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FMonomial.h"
#include "FixedAlgebra/FTps.h"
#include <complex>
#include <iosfwd>



// Template class FLieGenerator<T,N>.
// This code contains various tweaks for speed.
// ------------------------------------------------------------------------

template <class T, int N>
FLieGenerator<T, N>::
FLieGenerator(int order):
    itsOrder(order),
    bottomIndex(getBottomIndex(order)),
    topIndex(getTopIndex(order)),
    itsCoeffs(getSize(order)) {
    std::fill(itsCoeffs.begin(), itsCoeffs.end(), T(0));
}


template <class T, int N>
FLieGenerator<T, N>::
FLieGenerator(const FTps<T, 2 * N> &tps, int order):
    itsOrder(order),
    bottomIndex(getBottomIndex(order)),
    topIndex(getTopIndex(order)),
    itsCoeffs(getSize(order)) {
    if(order <= tps.getMaxOrder()) {
        std::copy(tps.begin() + bottomIndex, tps.begin() + topIndex,
                  itsCoeffs.begin());
    } else {
        std::fill(itsCoeffs.begin(), itsCoeffs.end(), T(0));
    }
}

template <class T, int N>
FLieGenerator<T, N>::
FLieGenerator():
    itsOrder(-1),
    bottomIndex(0),
    topIndex(0),
    itsCoeffs(1) {
    itsCoeffs[0] = T(0);
}


template <class T, int N>
FLieGenerator<T, N>::
FLieGenerator(const FLieGenerator &rhs):
    itsOrder(rhs.itsOrder),
    bottomIndex(rhs.bottomIndex),
    topIndex(rhs.topIndex),
    itsCoeffs(rhs.itsCoeffs)
{}


template <class T, int N>
FLieGenerator<T, N>::
~FLieGenerator()
{}


template <class T, int N>
const FLieGenerator<T, N> &FLieGenerator<T, N>::
operator=(const FLieGenerator &rhs) {
    itsCoeffs = rhs.itsCoeffs;
    itsOrder = rhs.itsOrder;
    bottomIndex = rhs.bottomIndex;
    topIndex = rhs.topIndex;
    return *this;
}


template<class T, int N> inline
T *FLieGenerator<T, N>::
begin() {
    return itsCoeffs.begin();
}


template<class T, int N> inline
const T *FLieGenerator<T, N>::
begin() const {
    return itsCoeffs.begin();
}


template<class T, int N> inline
T *FLieGenerator<T, N>::end() {
    return itsCoeffs.end();
}


template<class T, int N> inline
const T *FLieGenerator<T, N>::
end() const {
    return itsCoeffs.end();
}


template<class T, int N> inline
T &FLieGenerator<T, N>::operator[](int i) {
    return itsCoeffs[i-bottomIndex];
}


template<class T, int N> inline
const T &FLieGenerator<T, N>::operator[](int i) const {
    return itsCoeffs[i-bottomIndex];
}


template<class T, int N>
FLieGenerator<T, N> FLieGenerator<T, N>::operator-() const {
    if(itsOrder < 0) {
        return *this;
    } else {
        FLieGenerator<T, N> result(itsOrder);
        std::transform(itsCoeffs.begin(), itsCoeffs.end(),
                       result.itsCoeffs.begin(), std::negate<T>());
        return result;
    }
}


template<class T, int N>
FLieGenerator<T, N> &FLieGenerator<T, N>::operator*=(const T &val) {
    if(itsOrder < 0) {
        return *this;
    } else {
        std::transform(itsCoeffs.begin(), itsCoeffs.end(), itsCoeffs.begin(),
                       std::bind2nd(std::multiplies<T>(), val));
        return *this;
    }
}


template<class T, int N>
FLieGenerator<T, N> &FLieGenerator<T, N>::operator/=(const T &val) {
    if(itsOrder < 0) {
        return *this;
    } else {
        std::transform(itsCoeffs.begin(), itsCoeffs.end(), itsCoeffs.begin(),
                       std::bind2nd(std::divides<T>(), val));
        return *this;
    }
}


template<class T, int N>
FLieGenerator<T, N> &FLieGenerator<T, N>::
operator+=(const FLieGenerator<T, N> &rhs) {
    if(itsOrder == rhs.itsOrder) {
        std::transform(itsCoeffs.begin(), itsCoeffs.end(), rhs.itsCoeffs.begin(),
                       itsCoeffs.begin(), std::plus<T>());
    } else if(isZero()) {
        *this = rhs;
    } else if(! rhs.isZero()) {
        throw SizeError("operator+=(FLieGenerator)", "Inconsistent orders.");
    }

    return *this;
}


template<class T, int N>
FLieGenerator<T, N> &FLieGenerator<T, N>::
operator-=(const FLieGenerator<T, N> &rhs) {
    if(itsOrder == rhs.itsOrder) {
        std::transform(itsCoeffs.begin(), itsCoeffs.end(), rhs.itsCoeffs.begin(),
                       itsCoeffs.begin(), std::minus<T>());
    } else if(isZero()) {
        *this = rhs;
    } else if(! rhs.isZero()) {
        throw SizeError("operator+=(FLieGenerator)", "Inconsistent orders.");
    }

    return *this;
}


template <class T, int N>
void FLieGenerator<T, N>::
clear() {
    std::fill(itsCoeffs.begin(), itsCoeffs.end(), T(0));
}


template <class T, int N>
FLieGenerator<T, N> FLieGenerator<T, N>::derivative(int var) const {
    if(getOrder() > 0) {
        FLieGenerator<T, N> result(getOrder() - 1);
        const Array1D<int> &product = FTpsData<2 * N>::getProductArray(var + 1);

        for(int i = result.getBottomIndex(); i < result.getTopIndex(); ++i) {
            int k = product[i];
            result[i] = (*this)[k] * double(FTpsData<2 * N>::getExponents(k)[var]);
        }

        return result;
    } else {
        return FLieGenerator<T, N>();
    }
}


template <class T, int N>
bool FLieGenerator<T, N>::
isZero() const {
    if(itsOrder < 0) return true;

    for(const T *i = itsCoeffs.begin(); i != itsCoeffs.end(); ++i) {
        if(*i != T(0)) return false;
    }

    return true;
}


template <class T, int N>
FLieGenerator<T, N> FLieGenerator<T, N>::
scale(const FLieGenerator &rhs) const {
    FLieGenerator result(itsOrder);
    std::transform(itsCoeffs.begin(), itsCoeffs.end(), rhs.itsCoeffs.begin(),
                   result.itsCoeffs.begin(), std::multiplies<T>());
    return result;
}


template <class T, int N> template <class U>
FLieGenerator<U, N> FLieGenerator<T, N>::
transform(const FMatrix<U, 2 * N, 2 * N> &mat) const {
    int order = getOrder();
    FLieGenerator<U, N> trans(order);

    // Working array to accumulate partial products.
    Array1D<U> temp(getTopIndex());
    temp[0] = U(1);
    U *tt = temp.begin();

    const Array1D<TpsSubstitution> &table = FTpsData<2 * N>::getSubTable();
    for(int next = 1; next < table.size();) {
        const TpsSubstitution &s = table[next];
        int bot1 = getBottomIndex(s.order);
        int top1 = getTopIndex(s.order);
        for(int i = bot1; i < top1; ++i) tt[i] = U(0);
        const U *row = mat[s.variable];

        for(int v = 0; v < 2 * N; ++v) {
            U c = row[v];
            if(c != 0.0) {
                const Array1D<int> &prod = FTpsData<2 * N>::getProductArray(v + 1);
                int bot2 = getBottomIndex(s.order - 1);
                int top2 = getTopIndex(s.order - 1);
                for(int k = bot2; k < top2; ++k) {
                    tt[prod[k]] += c * tt[k];
                }
            }
        }

        if(s.order < order) {
            ++next;
        } else if(s.order == order) {
            T coeff = (*this)[s.index];
            if(coeff != T(0)) {
                for(int k = bot1; k < top1; ++k) {
                    trans[k] += coeff * tt[k];
                }
            }
            next = s.skip;
        } else {
            next = s.skip;
        }
    }

    return trans;
}


template <class T, int N>
inline int FLieGenerator<T, N>::
getOrder() const {
    return itsOrder;
}


template <class T, int N>
inline int FLieGenerator<T, N>::
getBottomIndex() const {
    return bottomIndex;
}


template <class T, int N>
inline int FLieGenerator<T, N>::
getTopIndex() const {
    return topIndex;
}


template <class T, int N>
inline int FLieGenerator<T, N>::
getSize(int order) {
    return FTpsData<2 * N>::getSize(order) - FTpsData<2 * N>::getSize(order - 1);
}


template <class T, int N>
inline int FLieGenerator<T, N>::
getBottomIndex(int order) {
    return FTpsData<2 * N>::getSize(order - 1);
}


template <class T, int N>
inline int FLieGenerator<T, N>::
getTopIndex(int order) {
    return FTpsData<2 * N>::getSize(order);
}


// Global functions.
// ------------------------------------------------------------------------

template <class T, int N>
FLieGenerator<T, N>
operator+(const FLieGenerator<T, N> &x, const FLieGenerator<T, N> &y) {
    FLieGenerator<T, N> z(x);
    return z += y;
}


template <class T, int N>
FLieGenerator<T, N>
operator-(const FLieGenerator<T, N> &x, const FLieGenerator<T, N> &y) {
    FLieGenerator<T, N> z(x);
    return z -= y;
}


template <class T, int N>
FLieGenerator<T, N>
operator*(const FLieGenerator<T, N> &x, const T &y) {
    FLieGenerator<T, N> z(x);
    return z *= y;
}


template <class T, int N>
FLieGenerator<T, N>
operator*(const T &x, const FLieGenerator<T, N> &y) {
    FLieGenerator<T, N> z(y);
    return z *= x;
}


template <class T, int N>
FLieGenerator<T, N>
operator*(const FLieGenerator<T, N> &x, const FLieGenerator<T, N> &y) {
    if(x.isZero()) {
        return FLieGenerator<T, N>();
    } else if(y.isZero()) {
        return FLieGenerator<T, N>();
    } else {
        FLieGenerator<T, N> z(x.getOrder() + y.getOrder());

        int bot1 = x.getBottomIndex();
        int top1 = x.getTopIndex();
        int bot2 = y.getBottomIndex();
        int top2 = y.getTopIndex();

        for(int i1 = bot1; i1 < top1; ++i1) {
            const Array1D<int> &p = FTpsData<2 * N>::getProductArray(i1);
            T f = x[i1];
            if(f != 0.0) {
                for(int i2 = bot2; i2 < top2; ++i2) {
                    z[p[i2]] += f * y[i2];
                }
            }
        }

        return z;
    }
}


template <class T, int N>
FLieGenerator<T, N>
operator/(const FLieGenerator<T, N> &x, const T &y) {
    FLieGenerator<T, N> z(x);
    return z /= y;
}


template <class T, int N>
FLieGenerator<T, N>
real(const FLieGenerator<std::complex<T>, N> &x) {
    FLieGenerator<T, N> z(x.getOrder());
    T *i1 = z.begin();
    const std::complex<T> *i2 = x.begin();

    while(i1 != z.end()) {
        *i1++ = std::real(*i2++);
    }

    return z;
}


template <class T, int N>
FLieGenerator<T, N>
imag(const FLieGenerator<std::complex<T>, N> &x) {
    FLieGenerator<T, N> z(x.getOrder());
    T *i1 = z.begin();
    const std::complex<T> *i2 = x.begin();

    while(i1 != z.end()) {
        *i1++ = std::imag(*i2++);
    }

    return z;
}


template <class T, int N>
FLieGenerator<std::complex<T>, N>
toComplex(const FLieGenerator<T, N> &x) {
    FLieGenerator<std::complex<T>, N> z(x.getOrder());
    std::complex<T> *i1 = z.begin();
    const T *i2 = x.begin();

    while(i1 != z.end()) {
        *i1++ = std::complex<T>(*i2++, T(0));
    }

    return z;
}


template <class T, int N> FLieGenerator<T, N>
PoissonBracket(const FLieGenerator<T, N> &f, const FLieGenerator<T, N> &g) {
    if(f.isZero()) {
        return FLieGenerator<T, N>();
    } else if(g.isZero()) {
        return FLieGenerator<T, N>();
    } else {
        // Index limits for the derivatives.
        int ord_f = f.getOrder() - 1;
        int ord_g = g.getOrder() - 1;
        int top_f = f.getBottomIndex();
        int top_g = g.getBottomIndex();
        int bot_f = (top_f * ord_f) / (2 * N + ord_f);
        int bot_g = (top_g * ord_g) / (2 * N + ord_g);
        const T *ff = f.begin() - f.getBottomIndex();
        const T *gg = g.begin() - g.getBottomIndex();

        FLieGenerator<T, N> h(ord_f + ord_g);
        T *hh = h.begin() - h.getBottomIndex();

        for(int ind_1 = bot_f; ind_1 < top_f; ++ind_1) {
            const Array1D<int> &prod_f = FTpsData<2 * N>::getProductArray(ind_1);
            const FMonomial<2 * N> &exp_f = FTpsData<2 * N>::getExponents(ind_1);

            for(int ind_g = bot_g; ind_g < top_g; ++ind_g) {
                const Array1D<int> &prod_g = FTpsData<2 * N>::getProductArray(ind_g);
                const FMonomial<2 * N> &exp_g = FTpsData<2 * N>::getExponents(ind_g);

                for(int var = 0; var < 2 * N; var += 2) {
                    // Index of product monomial.
                    int ind_h = prod_f[ind_g];

                    // The relevant coefficients of the four derivatives.
                    double f_q = ff[prod_f[var+1]] * T(exp_f[var]   + 1);
                    double g_p = gg[prod_g[var+2]] * T(exp_g[var+1] + 1);
                    double f_p = ff[prod_f[var+2]] * T(exp_f[var+1] + 1);
                    double g_q = gg[prod_g[var+1]] * T(exp_g[var]   + 1);

                    // Accumulate.
                    hh[ind_h] += (f_q * g_p - f_p * g_q);
                }
            }
        }

        return h;
    }
}


template <class T, int N>
std::ostream &operator<<(std::ostream &os, const FLieGenerator<T, N> &gen) {
    os << "Tps" << std::setw(4) << gen.getOrder() << std::setw(4)
       << gen.getOrder() << std::setw(4) << (2 * N) << std::endl;
    std::streamsize old_prec = os.precision(14);
    os.setf(std::ios::scientific, std::ios::floatfield);

    for(int i = gen.getBottomIndex(); i < gen.getTopIndex(); ++i) {
        if(gen[i] != T(0)) {
            os << std::setw(24) << gen[i];
            const FMonomial<2 * N> &index = FTpsData<2 * N>::getExponents(i);

            for(int var = 0; var < 2 * N; var++) {
                os << std::setw(3) << index[var];
            }

            os << std::endl;
        }
    }

    os << std::setw(24) << T(0);

    for(int var = 0; var < 2 * N; var++) {
        os << std::setw(3) << (-1);
    }

    os << std::endl;
    os.setf(std::ios::fixed, std::ios::floatfield);
    os.precision(old_prec);
    return os;
}

#endif
