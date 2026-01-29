#ifndef CLASSIC_FVps_CC
#define CLASSIC_FVps_CC

// ------------------------------------------------------------------------
// $RCSfile: FVps.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2.2.8 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template Class: FVps<T,N>
//   Vector of power series with dimension N and N variables.
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2004/11/18 22:18:06 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

//#define DEBUG_FVps_CC

#include "FixedAlgebra/FLUMatrix.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsData.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/LinearMap.h"
#include "FixedAlgebra/TransportMap.h"
#include "Utilities/ConvergenceError.h"
#include "Utilities/DomainError.h"
#include "Utilities/FormatError.h"
#include "Utilities/LogicalError.h"
#include "Utilities/CLRangeError.h"
#include <iostream>
#include <iterator>
#include <list>

// Template class FVps<T,N>
// ------------------------------------------------------------------------

template <class T, int N>
FVps<T, N>::FVps() {
    identity();
}


template <class T, int N>
FVps<T, N>::FVps(const FVps<T, N> &rhs) {
    for(int i = 0; i < N; ++i) data[i] = rhs.data[i];
}


template <class T, int N>
FVps<T, N>::FVps(const LinearMap<T, N> &rhs) {
    identity();
    for(int i = 0; i < N; ++i) {
        if(rhs[i][0] != T(0)) data[i].setMinOrder(0);
        for(int j = 0; j <= N; ++j) data[i][j] = rhs[i][j];
    }
}


template <class T, int N>
FVps<T, N>::FVps(const TransportMap<T, N> &rhs) {
    static const int SIZE = FTpsData<N>::getSize(2);
    for(int i = 0; i < N; ++i) {
        data[i] = FTps<T, N>(0, 2, 2);
        for(int j = 0; j < SIZE; ++j) data[i][j] = rhs[i][j];
    }
}


template <class T, int N>
FVps<T, N>::FVps(int minOrder, int maxOrder, int trcOrder) {
    for(int i = 0; i < N; ++i)
        data[i] = FTps<T, N>(minOrder, maxOrder, trcOrder);
}


template <class T, int N>
FVps<T, N>::FVps(const FMatrix<T, N, N> &x) {
    identity();
    for(int i = 0; i < N; ++i)
        for(int j = 0; j < N; j++) data[i][j+1] = x(i, j);
}


template <class T, int N>
FVps<T, N>::FVps(const FVector<T, N> &x) {
    for(int i = 0; i < N; i++) data[i] = FTps<T, N>(x[i]);
}


template <class T, int N>
FVps<T, N>::~FVps()
{}


template <class T, int N>
const FVps<T, N> &FVps<T, N>::operator=(const FVps<T, N> &rhs) {
    if(&rhs != this) for(int i = 0; i < N; ++i) data[i] = rhs.data[i];
    return *this;
}


template <class T, int N>
void FVps<T, N>::identity() {
    for(int i = 0; i < N; ++i) data[i] = FTps<T, N>::makeVariable(i);
}


template <class T, int N>
void FVps<T, N>::zero() {
    for(int i = 0; i < N; ++i) data[i] = T(0);
}


template <class T, int N>
const FTps<T, N> &FVps<T, N>::getComponent(int index) const {
    if(index < 0 || index >= N)
        throw CLRangeError("FVps::getComponent()", "Index out of range.");

    return data[index];
}


template <class T, int N>
void FVps<T, N>::setComponent(int index, const FTps<T, N> &value) {
    if(index < 0 || index >= N)
        throw CLRangeError("FVps::setComponent()", "Index out of range.");

    data[index] = value;
}


template <class T, int N> inline
const FTps<T, N> &FVps<T, N>::operator[](int index) const {
    return data[index];
}


template <class T, int N> inline
FTps<T, N> &FVps<T, N>::operator[](int index) {
    return data[index];
}


template <class T, int N>
int FVps<T, N>::getDimension() const {
    return N;
}


template <class T, int N>
int FVps<T, N>::getVariables() const {
    return N;
}


template <class T, int N>
int FVps<T, N>::getMinOrder() const {
    const FTps<T, N> *p = data + N - 1;
    int minOrder = p->getMinOrder();

    while(p-- > data) minOrder = std::min(minOrder, p->getMinOrder());
    return minOrder;
}


template <class T, int N>
void FVps<T, N>::setMinOrder(int order) {
    for(int i = 0; i < N; ++i) data[i].setMinOrder(order);
}


template <class T, int N>
int FVps<T, N>::getMaxOrder() const {
    const FTps<T, N> *p = data + N;
    int maxOrder = 0;

    while(p-- > data) maxOrder = std::max(maxOrder, p->getMaxOrder());
    return maxOrder;
}


template <class T, int N>
void FVps<T, N>::setMaxOrder(int order) {
    for(int i = 0; i < N; ++i) data[i].setMaxOrder(order);
}


template <class T, int N>
int FVps<T, N>::getTopOrder() const {
    const FTps<T, N> *p = data + N;
    int topOrder = 0;

    while(p-- > data) topOrder = std::max(topOrder, p->getMaxOrder());
    return topOrder;
}


template <class T, int N>
int FVps<T, N>::getTruncOrder() const {
    const FTps<T, N> *p = data + N - 1;
    int trcOrder = p->getTruncOrder();

    while(p-- > data) trcOrder = std::min(trcOrder, p->getTruncOrder());
    return trcOrder;
}


template <class T, int N>
void FVps<T, N>::setTruncOrder(int order) {
    for(int i = 0; i < N; ++i) data[i].setTruncOrder(order);
}


template <class T, int N>
FVps<T, N> FVps<T, N>::filter(int minOrder, int maxOrder, int trcOrder) const {
    // Default: trunc = FTps<T,N>::EXACT

    FVps<T, N> result;
    for(int i = 0; i < N; i++)
        result[i] = data[i].filter(minOrder, maxOrder, trcOrder);
    return result;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::truncate(int trunc) {
    return filter(0, trunc, trunc);
}


template <class T, int N>
FVps<T, N> FVps<T, N>::operator+() const {
    return *this;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::operator-() const {
    FVps<T, N> result;
    for(int i = 0; i < N; i++) result[i] = - data[i];
    return result;
}


template <class T, int N>
FVps<T, N>& FVps<T, N>::operator+=(const FVps<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] += rhs[i];
    return *this;
}


template <class T, int N>
FVps<T, N>& FVps<T, N>::operator-=(const FVps<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] -= rhs[i];
    return *this;
}


template <class T, int N>
FVps<T, N>& FVps<T, N>::operator+=(const FVector<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] += rhs[i];
    return *this;
}


template <class T, int N>
FVps<T, N>& FVps<T, N>::operator-=(const FVector<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] -= rhs[i];
    return *this;
}


template <class T, int N>
FVps<T, N> &FVps<T, N>::operator*=(const FTps<T, N> &rhs) {
    for(int i = 0; i < N; i++) data[i] *= rhs;
    return *this;
}

template <class T, int N>
FVps<T, N> FVps<T, N>::operator*(const FVps<T, N>& rhs) const {
    FVps<T, N> result;
    // go through variables (N) and compute truncated power series for them
    for(int i = 0; i < N; ++i) {
        // get truncated power series for variable i
        FTps<T, N> tps = this->getComponent(i);
        
        // initialize tps of result
        FTps<T, N> r = 0.0;
        // fake order --> gets updated inside this function
        r.setMinOrder(1);
        
        /* get coefficients and exponents of their monomials and multiply
         * truncated power series of rhs with appropriate power and coefficient
         */
        std::list<int> coeffs = tps.getListOfNonzeroCoefficients();
        
        for(std::list<int>::iterator it = coeffs.begin(); it != coeffs.end(); ++it) {
            
            FArray1D<int, N> expons = tps.extractExponents(*it);
            
            // represents the monomial --> is polynomial due to multiplication of each variable's polynomial
            FTps<T, N> mono = 1.0;
            
            for(int j = 0; j < N; ++j) {
                
                if (expons[j] != 0) {
                    // multiply each variable of the monomial of appropriate power, i.e. build monomial
                    FTps<T, N> tmp = rhs.getComponent(j);
                    mono = mono.multiply(tmp.makePower(expons[j]), FTps<T, N>::getGlobalTruncOrder());
                }
            }
            // multiply truncated power series with appropriate coefficient
            mono *= tps.getCoefficient(*it);
            
            /* sum up all polynomials that build the FTps of the result map for that variable,
             * make sure that the minimum order is correct.
             */
            r.setMinOrder(std::min(r.getMinOrder(), mono.getMinOrder()));
            r += mono;
        }        
        // computation of Tps of variable i finished
        result.setComponent(i,r);
    }
    
    return result;
}


template <class T, int N>
FVps<T, N> &FVps<T, N>::operator/=(const FTps<T, N> &rhs) {
    FTps<T, N> t = rhs.inverse();
    for(int i = 0; i < N; i++) data[i] *= t;
    return *this;
}


template <class T, int N>
FVps<T, N> &FVps<T, N>::operator*=(const T &rhs) {
    for(int i = 0; i < N; i++) data[i] *= rhs;
    return *this;
}


template <class T, int N>
FVps<T, N> &FVps<T, N>::operator/=(const T &rhs) {
    for(int i = 0; i < N; i++) data[i] /= rhs;
    return *this;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::inverse(int trunc) const {
    // Default: trunc = FTps<T,N>::EXACT

    // Get orders.
    int minOrder = getMinOrder(), maxOrder = getMaxOrder(), trcOrder = getTruncOrder();
    maxOrder = std::min(maxOrder, trunc);
    trcOrder = std::min(trcOrder, trunc);

    // Check sanity.
    if(maxOrder < minOrder) {
        std::cerr << " <*** ERROR ***> in FVps::inverse():\n";
        throw LogicalError("FVps<T,N>::inverse()", "Map truncated to a zero map.");
    }

    // Exceptions for non-invertible maps.
    if(minOrder > 1) {
        std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                  << "    Cannot invert a purely nonlinear map." << std::endl;
        throw DomainError("FVps<T,N>::inverse()");
    } else if(minOrder == 1) {
        if(maxOrder > 1 && trcOrder == FTps<T, N>::EXACT) {
            std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                      << "    Cannot invert an EXACT nonlinear map." << std::endl;
            throw DomainError("FVps<T,N>::inverse()");
        }
    } else { // minOrder == 0
        if(maxOrder == 0) {
            std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                      << "    Cannot invert a constant map." << std::endl;
            throw DomainError("FVps<T,N>::inverse()");
        }
        if(maxOrder > 1) {
            std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                      << "    Cannot invert a nonlinear map containing a constant term." << std::endl;
            throw DomainError("FVps<T,N>::inverse()");
        }
        if(trcOrder != FTps<T, N>::EXACT) {
            std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                      << "    Cannot invert a map with both constant and linear terms unless it is EXACT."
                      << std::endl;
            throw DomainError("FVps<T,N>::inverse()");
        }
    }

    // Invert linear part.
    FMatrix<T, N, N> t1inv;
    FVps<T, N> r1;
    try {
        FLUMatrix<T, N> lu(linearTerms());
        t1inv = lu.inverse();
        r1 = FVps<T, N>(t1inv);
    } catch(SingularMatrixError &smx) {
        std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                  << "    Cannot invert a map having a singular linear part." << std::endl;
        throw DomainError("FVps<T,N>::inverse()");
    }

    // Return linear case.
    if(maxOrder == 1) {
        if(minOrder == 0) r1 -= t1inv * constantTerm();
        r1.setTruncOrder(trcOrder);
        return r1;
    }

    // General case (minOrder == 1, maxOrder > 1, trcOrder != EXACT).
    // Inverse map computed order by order using the relations
    //   intitial R = R_1 = T_1^{-1};
    //   iterate  R = R_1 o (I - T_{2..m} o R)  trcOrder
    FVps<T, N> id;
    FVps<T, N> result = r1;
    FVps<T, N> T2n = filter(2, maxOrder, trcOrder);
    for(int m = 2; m <= trcOrder; ++m) {
        FVps<T, N> tr = T2n.substitute(result, m);
        result = t1inv * (id - tr);
    }

    return result;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::myInverse(int trunc) const {
    // Default: trunc = FTps<T,N>::EXACT

    // Get orders.
    int minOrder = getMinOrder(), maxOrder = getMaxOrder(), trcOrder = getTruncOrder();
    maxOrder = std::min(maxOrder, trunc);
    trcOrder = std::min(trcOrder, trunc);

    // Check sanity.
    if(maxOrder < minOrder) {
        std::cerr << " <*** ERROR ***> in FVps::inverse():\n";
        throw LogicalError("FVps<T,N>::inverse()", "Map truncated to a zero map.");
    }

    // Exceptions for non-invertible maps.
    if(minOrder > 1) {
        std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                  << "    Cannot invert a purely nonlinear map." << std::endl;
        throw DomainError("FVps<T,N>::inverse()");
    } else if(minOrder == 1) {
        if(maxOrder > 1 && trcOrder == FTps<T, N>::EXACT) {
            std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                      << "    Cannot invert an EXACT nonlinear map." << std::endl;
            throw DomainError("FVps<T,N>::inverse()");
        }
    } else { // minOrder == 0
        if(maxOrder == 0) {
            std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                      << "    Cannot invert a constant map." << std::endl;
            throw DomainError("FVps<T,N>::inverse()");
        }
        if(maxOrder > 1) {
            std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                      << "    Cannot invert a nonlinear map containing a constant term." << std::endl;
            throw DomainError("FVps<T,N>::inverse()");
        }
        if(trcOrder != FTps<T, N>::EXACT) {
            std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                      << "    Cannot invert a map with both constant and linear terms unless it is EXACT."
                      << std::endl;
            throw DomainError("FVps<T,N>::inverse()");
        }
    }

    // Invert linear part.
    FMatrix<T, N, N> t1inv;
    FVps<T, N> r1;
    try {
        FLUMatrix<T, N> lu(linearTerms());
        t1inv = lu.inverse();
        r1 = FVps<T, N>(t1inv);
    } catch(SingularMatrixError &smx) {
        std::cerr << " <*** ERROR ***> in FVps::inverse():\n"
                  << "    Cannot invert a map having a singular linear part." << std::endl;
        throw DomainError("FVps<T,N>::inverse()");
    }

    // Return linear case.
    if(maxOrder == 1) {
        if(minOrder == 0) r1 -= t1inv * constantTerm();
        r1.setTruncOrder(trcOrder);
        return r1;
    }

    // General case (minOrder == 1, maxOrder > 1, trcOrder != EXACT).
    // Inverse map computed order by order using the relations
    //   intitial R = R_1 = T_1^{-1};
    //   iterate  R = R_1 o (I - T_{2..m} o R)  trcOrder
    FVps<T, N> id;
    FVps<T, N> result = r1;
    FVps<T, N> T2n = filter(2, maxOrder, trcOrder);
    for(int m = 2; m <= trcOrder; ++m) {
        FVps<T, N> tr = T2n.substitute(result, m);
        result = t1inv * (id - tr);
        tr = substitute(result, m);
        result += t1inv * (id - tr);
    }

    return result;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::derivative(int var) const {
    FVps<T, N> result;
    for(int i = 0; i < N; i++) result[i] = data[i].derivative(var);
    return result;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::integral(int var) const {
    FVps<T, N> result;
    for(int i = 0; i < N; i++) result[i] = data[i].integral(var);
    return result;
}


template <class T, int N>
FVector<T, N> FVps<T, N>::constantTerm() const {
    FVector<T, N> result;
    for(int i = 0; i < N; i++) {
        if(data[i].getMinOrder() == 0) result[i] = data[i][0];
        else result[i] = T(0);
    }
    return result;
}


template <class T, int N>
FVector<T, N> FVps<T, N>::constantTerm(const FVector<T, N> &P) const {
    // Evaluate each component.
    FVector<T, N> result;
    for(int v = N; v-- > 0;)
        result[v] = (*this)[v].evaluate(P);
    return result;
}


template <class T, int N>
FMatrix<T, N, N> FVps<T, N>::linearTerms() const {
    FMatrix<T, N, N> result;
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++) result(i, j) = data[i][j+1];
    return result;
}


template <class T, int N>
FMatrix<T, N, N> FVps<T, N>::linearTerms(const FVector<T, N> &P) const {
    // Evaluate monomials.
    int maxOrder = getMaxOrder();
    if(maxOrder) --maxOrder;
    Array1D<T> monoms = FTps<T, N>::evalMonoms(P, maxOrder);
    T *m = monoms.begin();

    FMatrix<T, N, N> result;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            FTps<T, N> gzij = data[i].derivative(j);
            int ks = FTps<T, N>::orderStart(gzij.getMinOrder());
            int ke = FTps<T, N>::orderEnd(gzij.getMaxOrder());
            T *dk = gzij.begin() + ks, *mk = m + ks, *mke = m + ke;
            T rij = T(0);
            while(mk != mke) rij += *dk++ * *mk++;
            result(i, j) = rij;
        }
    }
    return result;
}


template <class T, int N>
Array1D<int> FVps<T, N>::getSubstOrders(const FVps<T, N> &rhs, int trunc) const {
    // Defaul: trunc = EXACT

    // Get orders.
    Array1D<int> ordersL(3), ordersR(3);
    ordersL[0] = getMinOrder(),     ordersL[1] = getMaxOrder(),     ordersL[2] = getTruncOrder();
    ordersR[0] = rhs.getMinOrder(), ordersR[1] = rhs.getMaxOrder(), ordersR[2] = rhs.getTruncOrder();

    Array1D<int> result = FTps<T, N>::getSubstOrders(ordersL, ordersR, trunc);
    return result;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::substitute(const FMatrix<T, N, N> &mat, int n) const {
    // Check sanity.
    if(n < 0)
        throw LogicalError("FVps<T,N>::substitute(mat,n)",
                           "Transformation order, n, is negative.");
    else if(n > FTps<T, N>::getGlobalTruncOrder())
        throw LogicalError("FVps<T,N>::substitute(mat,n)",
                           "Transformation order, n, exceeds globalTruncOrder.");

    // Get orders; if necessary, make LHS have uniform min, max, and trc orders.
    FVps<T, N> f = *this;
    int minOrder = getMinOrder(), maxOrder = getMaxOrder(), trcOrder = getTruncOrder();
    for(int k = N; k-- > 0;) {
        if(f[k].getMinOrder()   != minOrder) f[k].setMinOrder(minOrder);
        if(f[k].getMaxOrder()   != maxOrder) f[k].setMaxOrder(maxOrder);
        if(f[k].getTruncOrder() != trcOrder) f[k].setTruncOrder(trcOrder);
    }

    //Allocate result.
    FVps<T, N> result(minOrder, maxOrder, trcOrder);
    for(int k = N; k-- > 0;)
        std::copy(f[k].begin(minOrder), f[k].end(maxOrder), result[k].begin(minOrder));

    // Return trivial cases.
    if(n > trcOrder) {
#ifdef DEBUG_FVps_CC
        std::cerr << " <*** WARNING ***> from FTps<T,N>::substitute(mat,n):\n"
                  << "    Transformation order exceeds truncation order;\n"
                  << "    returning map unchanged." << std::endl;
#endif
        return result;
    }
    if(n == 0 || n < minOrder || maxOrder < n) return result;

    // Allocate working array; use static
    // local memory to avoid fragmentation.
    static T *t = 0;
    static int max_n = -1;
    if(n > max_n) {
        if(t) delete [] t;
        t = new T[FTps<T, N>::getSize(n)];
        max_n = n;
    }

    // Initialisations.
    T *t1 = t + 1;
    const T *fj[N];
    T *g[N];
    const Array1D<int> *oldvrbl = 0;
    int start_n = FTps<T, N>::orderStart(n), end_n = FTps<T, N>::orderEnd(n);
    for(int k = N; k-- > 0;) {
        fj[k] = f[k].begin(n);
        g[k] = result[k].begin();
        std::fill(g[k] + start_n, g[k] + end_n, T(0));
    }

    // Loop over order n monomials.
    for(int j = start_n; j < end_n; ++j) {
        // Skip monomials with coefficient zero.
        bool zeroQ = true;
        for(int k = N; k-- > 0;)
            if(*fj[k] != T(0)) zeroQ = false;
        if(zeroQ) {
            for(int k = N; k-- > 0;) ++fj[k];
            continue;
        }

        // Get current monomial's variable list; compare with old variable list.
        const Array1D<int> *vrbl = &FTpsData<N>::getVariableList(j);
        int vi = 0;
        if(oldvrbl)
            while((*vrbl)[vi] == (*oldvrbl)[vi]) ++vi;

        const T *mv;
        int ord;
        // If vi = 0, we must start at the beginning; otherwise,
        // we may re-use the first vi orders stored in t.
        if(vi == 0) {
            mv = mat[(*vrbl)[0]];
            std::copy(mv, mv + N, t1);
            ord = 2;
        } else ord = vi + 1;

        // In working array t, clear orders we can't use.
        std::fill(t + FTps<T, N>::orderStart(ord), t + end_n, T(0));
        // Build the remainder.
        while(ord <= n) {
            // Build next order part of transformed monomial by multiplying
            // the part that is one order lower by the transformed version
            // of the next variable in the variable list.
            int ord1 = ord - 1;
            int start_l = FTps<T, N>::orderStart(ord1), end_l = FTps<T, N>::orderEnd(ord1);
            mv = mat[(*vrbl)[ord1]];  // transformed version of next variable
            for(int k = 0; k < N; k++) {
                T mvk = mv[k];
                if(mvk == T(0)) continue;
                const Array1D<int> &prod = FTpsData<N>::getProductArray(k + 1);
                for(int l = start_l; l < end_l; l++) t[prod[l]] += mvk * t[l];
            }
            ++ord;
        }
        //Increment g[k] by fj[k] * transformed monomial.
        for(int k = N; k-- > 0;) {
            T *gk = g[k];
            T fjk = *fj[k];
            if(fjk != T(0))
                for(int i = start_n; i < end_n; i++) gk[i] += fjk * t[i];
        }
        // Save variable list for comparison with the next one.
        oldvrbl = vrbl;

        // Increment array of monomial pointers.
        for(int k = N; k-- > 0;) ++fj[k];
    }

    return result;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::substitute(const FMatrix<T, N, N> &mat, int nl, int nh) const {
    // Check sanity.
    if(nl > nh)
        throw LogicalError("FVps<T,N>::substitute(mat,nl,nh)",
                           "Inconsistent transformation orders: nl > nh.");
    if(nl < 0)
        throw LogicalError("FVps<T,N>::substitute(mat,nl,nh)",
                           "Transformation order nl is negative.");
    else if(nh > FTps<T, N>::getGlobalTruncOrder())
        throw LogicalError("FVps<T,N>::substitute(mat,nl,nh)",
                           "Transformation order nh exceeds globalTruncOrder.");

    // Get orders; if necessary, make LHS have uniform min, max, and trc orders..
    FVps<T, N> f = *this;
    int minOrder = getMinOrder(), maxOrder = getMaxOrder(), trcOrder = getTruncOrder();
    for(int k = N; k-- > 0;) {
        if(f[k].getMinOrder()   != minOrder) f[k].setMinOrder(minOrder);
        if(f[k].getMaxOrder()   != maxOrder) f[k].setMaxOrder(maxOrder);
        if(f[k].getTruncOrder() != trcOrder) f[k].setTruncOrder(trcOrder);
    }

    //Allocate result.
    FVps<T, N> result(minOrder, maxOrder, trcOrder);
    for(int k = N; k-- > 0;)
        std::copy(f[k].begin(minOrder), f[k].end(maxOrder), result[k].begin(minOrder));

    if(nh > trcOrder) {
#ifdef DEBUG_FVps_CC
        std::cerr << " <*** WARNING ***> from FVps<T,N>::substitute(mat,nl,nh):\n"
                  << "    Transformation order nh exceeds truncation order;\n"
                  << "    truncation order unchanged." << std::endl;
#endif
    }

    // Return trivial cases.
    if(nh == 0 || nh < minOrder || maxOrder < nl) return result;

    // Set and clear actual range of orders to transform.
    if(nl == 0) nl = 1;
    nl = std::max(nl, minOrder);
    nh = std::min(nh, maxOrder);
    for(int k = N; k-- > 0;)
        std::fill(result[k].begin(nl), result[k].end(nh), T(0));

    // Allocate working arrays; use static
    // local memory to avoid fragmentation.
    static T *t = 0;
    static int max_nh = -1;
    if(nh > max_nh) {
        if(t) delete [] t;
        t = new T[FTps<T, N>::getSize(nh)];
        max_nh = nh;
    }
    T *t1 = t + 1;

    // Initialisations.
    // Array element fp[k][m] points to the next order m monomial
    // to transform in k-th component.
    const T *fp[N][nh+1];
    T *g[N];
    for(int k = N; k-- > 0;) {
        for(int m = nl; m <= nh; ++m) fp[k][m] = f[k].begin(m);
        g[k] = result[k].begin();
    }
    const Array1D<int> *oldvrbl = 0;
    int start_nh = FTps<T, N>::orderStart(nh), end_nh = FTps<T, N>::orderEnd(nh);
    int nh1 = nh - 1, nh2 = nh - 2;

    // Loop over order nh monomials; construct lower orders along the way.
    for(int j = start_nh; j < end_nh; ++j) {
        // Get current monomial's variable list; compare with old variable list.
        const Array1D<int> *vrbl = &FTpsData<N>::getVariableList(j);
        int vk = 0;
        if(oldvrbl)
            while((*vrbl)[vk] == (*oldvrbl)[vk]) ++vk;

        // Determine which monomial pointers we shall need to increment.
        int jl = (*vrbl)[nh1], ni = nh2;
        while(ni >= 0 && (*vrbl)[ni] == jl) --ni;
        ni += 2;
        ni = std::max(ni, nl);
        // Determine which monomials contribute this round.
        int n1 = std::max(nl, ni), n2 = nh;
        while(n1 <= n2) {
            bool zeroQ = true;
            for(int k = N; k-- > 0;) {
                if(*fp[k][n1] != T(0)) {
                    zeroQ = false;
                    break;
                }
            }
            if(zeroQ) ++n1;
            else break;
        }
        while(n2 > n1) {
            bool zeroQ = true;
            for(int k = N; k-- > 0;) {
                if(*fp[k][n2] != T(0)) {
                    zeroQ = false;
                    break;
                }
            }
            if(zeroQ) --n2;
            else break;
        }
        // Skip if all monomials have coefficient zero.
        if(n1 > n2) {
            for(int k = N; k-- > 0;)
                for(int m = ni; m <= nh; ++m) ++fp[k][m];
            continue;
        }

        const T *mv;
        int ord;
        // If vk = 0, we must start at the beginning; otherwise,
        // we may keep the first vk orders stored in t.
        if(vk == 0) {
            mv = mat[(*vrbl)[0]];
            std::copy(mv, mv + N, t1);
            ord = 2;
        } else ord = vk + 1;

        // In working array t, clear orders we can't use.
        std::fill(t + FTps<T, N>::orderStart(ord), t + end_nh, T(0));

        // Build the remainder.
        while(ord <= nh) {
            // Build next order part of transformed monomial by multiplying
            // the part that is one order lower by the transformed version
            // of the next variable in the variable list.
            int ord1 = ord - 1;
            int start_l = FTps<T, N>::orderStart(ord1), end_l = FTps<T, N>::orderEnd(ord1);
            mv = mat[(*vrbl)[ord1]];  // transformed version of next variable
            for(int k = 0; k < N; k++) {
                T mvk = mv[k];
                if(mvk == T(0)) continue;
                const Array1D<int> &prod = FTpsData<N>::getProductArray(k + 1);
                for(int l = start_l; l < end_l; ++l) t[prod[l]] += mvk * t[l];
            }
            ++ord;
        }
        // Increment g[k] by f[k][j] * transformed monomial.
        // and increment pointers in fp[][].
        for(int k = N; k-- > 0;) {
            for(int m = n1; m <= n2; ++m) {
                const T fkj = *fp[k][m];
                int start_m = FTps<T, N>::orderStart(m), end_m = FTps<T, N>::orderEnd(m);
                for(int i = start_m; i < end_m; i++) g[k][i] += fkj * t[i];
            }
            for(int m = ni; m <= nh; ++m) ++fp[k][m];
        }

        // Save variable list for comparison with the next one.
        oldvrbl = vrbl;
    }

    return result;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::substitute(const FMatrix<T, N, N> &mat) const {
    return substitute(mat, getMinOrder(), getMaxOrder());
}


template <class T, int N>
FVps<T, N> FVps<T, N>::substitute(const FVps<T, N> &rhs, int trunc) const {
    // Default: trunc = FTps<T,N>::EXACT

    // Get orders; if necessary, make LHS have uniform min, max, and trc orders..
    FVps<T, N> f = *this;
    int f_min = getMinOrder(), f_max = getMaxOrder(), f_trc = getTruncOrder();
    Array1D<int> orders = getSubstOrders(rhs, trunc);
    int g_min = orders[0], g_max = orders[1], g_trc = orders[2];
    for(int k = N; k-- > 0;) {
        if(f[k].getMinOrder()   != f_min) f[k].setMinOrder(f_min);
        if(f[k].getMaxOrder()   != f_max) f[k].setMaxOrder(f_max);
        if(f[k].getTruncOrder() != f_trc) f[k].setTruncOrder(f_trc);
    }

    // Make sure we don't trip over globalTruncOrder.
    if(g_trc != FTps<T, N>::EXACT && g_trc > FTps<T, N>::getGlobalTruncOrder())
        throw LogicalError("FVps::substitute(FVps rhs, int trunc)",
                           "Truncation order exceeds globalTruncOrder!");

    // Return trivial case.
    if(g_min > g_max) return FVps<T, N>(g_trc, g_trc, g_trc);

    //Allocate result.
    FVps<T, N> result(g_min, g_max, g_trc);
    if(f_min == 0)
        for(int k = N; k-- > 0;) result[k][0] = *f[k].begin();
    if(f_max == 0) return result;

    // Set actual range of orders to transform
    int nl = f_min, nh = f_max;
    if(nl == 0) nl = 1;

    // Allocate working arrays.
    const T *fp[N][nh+1];
    Array1D< FTps<T, N> > t(nh + 1);

    // Initialisations.
    // Array element fp[k][m] points to the next order m monomial
    // to transform in the k-th component.
    for(int k = N; k-- > 0;)
        for(int m = nl; m <= nh; ++m) fp[k][m] = f[k].begin(m);
    const Array1D<int> *oldvrbl = 0;
    int start_nh = FTps<T, N>::orderStart(nh), end_nh = FTps<T, N>::orderEnd(nh);
    int nh1 = nh - 1, nh2 = nh - 2;

    // Loop over order nh monomials; construct lower orders along the way.
    for(int j = start_nh; j < end_nh; ++j) {
        // Get current monomial's variable list; compare with old variable list.
        const Array1D<int> *vrbl = &FTpsData<N>::getVariableList(j);
        int vk = 0;
        if(oldvrbl)
            while((*vrbl)[vk] == (*oldvrbl)[vk]) ++vk;

        // Determine which monomial pointers we shall need to increment.
        int jl = (*vrbl)[nh1], ni = nh2;
        while(ni >= 0 && (*vrbl)[ni] == jl) --ni;
        ni += 2;
        ni = std::max(ni, nl);
        // Determine which monomials contribute this round.
        int n1 = std::max(nl, ni), n2 = nh;
        while(n1 <= n2) {
            bool zeroQ = true;
            for(int k = N; k-- > 0;) {
                if(*fp[k][n1] != T(0)) {
                    zeroQ = false;
                    break;
                }
            }
            if(zeroQ) ++n1;
            else break;
        }
        while(n2 > n1) {
            bool zeroQ = true;
            for(int k = N; k-- > 0;) {
                if(*fp[k][n2] != T(0)) {
                    zeroQ = false;
                    break;
                }
            }
            if(zeroQ) --n2;
            else break;
        }
        // Skip if all monomials have coefficient zero.
        if(n1 > n2) {
            for(int k = N; k-- > 0;)
                for(int m = ni; m <= nh; ++m) ++fp[k][m];
            continue;
        }

        // If vk = 0, we must start at the beginning; otherwise,
        // we may keep the first vk orders stored in t.
        int ord;
        if(vk == 0) {
            t[1] = rhs[(*vrbl)[0]];
            ord = 2;
        } else ord = vk + 1;

        // Build the remainder.
        while(ord <= nh) {
            // Build next order part of transformed monomial by multiplying
            // the part that is one order lower by the transformed version
            // of the next variable in the variable list.
            int ord1 = ord - 1;
            t[ord] = t[ord1].multiply(rhs[(*vrbl)[ord1]], g_trc);
            ++ord;
        }

        // Increment result by f[k][j] * transformed monomial,
        // and increment pointers in fp[].
        for(int k = N; k-- > 0;) {
            const T **fpk = fp[k];
            for(int m = n1; m <= n2; ++m) result[k] += *fpk[m] * t[m];
            for(int m = ni; m <= nh; ++m) ++fpk[m];
            //for (int m = n1; m <= n2; ++m) result[k] += *fpk[m] * t[m];
            //for (int m = ni; m <= nh; ++m) ++fpk[m];
        }

        // Save variable list for comparison with the next one.
        oldvrbl = vrbl;
    }

    return result;
}


template <class T, int N>
FVps<T, N> FVps<T, N>::substituteInto(const FMatrix<T, N, N> &lhs) const {
    FVps<T, N> result;

    for(int i = 0; i < N; ++i) {
        FTps<T, N> sum = lhs(i, 0) * data[0];
        for(int j = 1; j < N; ++j) sum += lhs(i, j) * data[j];
        result[i] = sum;
    }

    return result;
}

template <class T, int N>
FTps<T, N> FVps<T, N>::getFTps(const FArray1D<int, N>& power) const {
    
    // function does not handle negative powers   
    if ( std::any_of(power.begin(), power.end(), [&](int p) { return p < 0; }) )
        throw LogicalError("FVps<T,N>::getFTps(power)", "Negative power.");
    
    // initial Tps
    FTps<T, N> result = 1.0;
    
    // go through variables and multiply its power to "result"
    for (int i = 0; i < N; ++i) {
        // get polynomial
        FTps<T, N> rhs = getComponent(i);
        
        // multiply polynomials
        for (int j = 0; j < power[i]; ++j)
            result = result.multiply(rhs, FTps<T, N>::getGlobalTruncOrder()); // make sure that no global trunc exceeding
    }
    
    return result;/*.truncate(FTps<T, N>::getGlobalTruncOrder());*/
}

template <class T, int N>
std::istream &FVps<T, N>::get(std::istream &is) {
    is.flags(std::ios::skipws);
    char head[4];
    (is >> std::ws).get(head, 4);

    if(strcmp(head, "FVps") != 0)
        throw FormatError("FVps::get()", "Flag word \"FVps\" missing.");

    int nDim;
    is >> nDim;
    if(nDim != N) throw FormatError("FVps::get()", "Invalid FVps dimension");

    // Read into temporary for exception safety.
    FVps<T, N> result;
    for(int i = 0; i < N; i++) is >> result.data[i];
    *this = result;
    return is;
}


template <class T, int N>
std::ostream &FVps<T, N>::put(std::ostream &os) const {
    os << "FVps " << N << std::endl;
    for(int i = 0; i < N; i++) os << data[i];
    return os;
}


// Global Operators on FVps<T,N>
// ------------------------------------------------------------------------

template <class T, int N>
FVps<T, N> operator+(const FVps<T, N> &lhs, const FVps<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] + rhs[i];
    return result;
}


template <class T, int N>
FVps<T, N> operator-(const FVps<T, N> &lhs, const FVps<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] - rhs[i];
    return result;
}


template <class T, int N>
FVps<T, N> operator+(const FVps<T, N> &lhs, const FVector<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] + rhs[i];
    return result;
}


template <class T, int N>
FVps<T, N> operator-(const FVps<T, N> &lhs, const FVector<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] - rhs[i];
    return result;
}


template <class T, int N>
FVps<T, N> operator+(const FVector<T, N> &lhs, const FVps<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] + rhs[i];
    return result;
}


template <class T, int N>
FVps<T, N> operator-(const FVector<T, N> &lhs, const FVps<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] - rhs[i];
    return result;
}


template <class T, int N>
FVector<T, N> operator*(const FVps<T, N> &lhs, const FVector<T, N> &rhs) {
    FVector<T, N> result;
    for (int i = 0; i < N; ++i) {
        FTps<T, N> tps = lhs.getComponent(i);
        result[i] = tps.evaluate(rhs);
    }
    return result;
}


template <class T, int N>
FVps<T, N> operator*(const FVps<T, N> &lhs, const FTps<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] * rhs;
    return result;
}


template <class T, int N>
FVps<T, N> operator*(const FTps<T, N> &lhs, const FVps<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs * rhs[i];
    return result;
}


template <class T, int N>
FVps<T, N> operator*(const FVps<T, N> &lhs, const T &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] * rhs;
    return result;
}


template <class T, int N>
FVps<T, N> operator*(const T &lhs, const FVps<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs * rhs[i];
    return result;
}

template <class T, int N>
FVps<T, N> operator*(const FMatrix<T, N, N> &lhs, const FVps<T, N> &rhs) {
    return rhs.substituteInto(lhs);
}


template <class T, int N>
FVps<T, N> operator/(const FVps<T, N> &lhs, const FTps<T, N> &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] / rhs;
    return result;
}


template <class T, int N>
FVps<T, N> operator/(const FVps<T, N> &lhs, const T &rhs) {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = lhs[i] / rhs;
    return result;
}


template <class T, int N> FVps<T, N>
ExpMap(const FTps<T, N> &H, const FVps<T, N> &map, int trunc) {
    //std::cerr << "==> In ExpMap(H,map,trunc)" << std::endl;
    // Default: trunc = FTps<T,N>::EXACT

    // Limit number of iterations.
    const int MAX_ITER = 400;

    // We really ought to throw an exception if H contains linear terms and is not exact,
    // but we're just going to complain!!
    bool FD = false;
    if(H.getTruncOrder() != FTps<T, N>::EXACT && H.getMinOrder() < 2) {
        FD = true;
#ifdef DEBUG_FVps_CC
        std::cerr << " <*** WARNING ***> from ExpMap(H,map,trunc):\n"
                  << "    Incomplete computation of feed-down terms.\n" << std::endl;
#endif
    }
    int fd_trc = std::min(H.getTruncOrder() - 1, map.getTruncOrder());

    // Construct dH = grad(H).J, s.t. :H:f = dH.grad(f).
    FVps<T, N> dH;
    for(int i = 0; i < N; i += 2) {
        dH[i]   = - H.derivative(i + 1);
        dH[i+1] =   H.derivative(i);
    }

    // Allocate result.
    FVps<T, N> expHmap;

    // Apply exp(:H:) to each component of map.
    for(int var = 0; var < N; var++) {
        // Initialize variables.
        FTps<T, N> expHf = map[var];
        FTps<T, N> dHkf = map[var];
        FTps<T, N> old = T(0);
        // Compute series; quit loop if we added nothing last time through.
        for(int k = 1; expHf != old; ++k) {
            if(k > MAX_ITER) {
                std::cerr << " present error:\n" << expHf - old << std::endl;
                throw ConvergenceError("ExpMap(const FTps<T,N> &H, const FVps<T,N> &map)",
                                       "No convergence in ExpMap(H,map)");
            }
            // Don't initialize dHk1f to 0, as that sets minOrder to 0!
            old = expHf;
            FTps<T, N> ddHkf = dHkf.derivative(0);
            if(FD) ddHkf.setTruncOrder(fd_trc);
            FTps<T, N> dHk1f = dH[0].multiply(ddHkf, trunc);
            for(int v = 1; v < N; ++v) {
                FTps<T, N> ddHkf = dHkf.derivative(v);
                if(FD) ddHkf.setTruncOrder(fd_trc);
                dHk1f += dH[v].multiply(ddHkf, trunc);
                //dHk1f += dH[v].multiply(dHkf.derivative(v),trunc);
            }
            dHkf = dHk1f / T(k); // :H:^{k}/(k!)
            expHf += dHkf;
        }
        expHmap[var] = expHf;
    }

    //std::cerr << "==> Leaving ExpMap(H,map,trunc)" << std::endl;
    return expHmap;
}


template <class T, int N> FVps<T, N>
PoissonBracket(const FTps<T, N> &x, const FVps<T, N> &y, int trunc) {
    FVps<T, N> z;
    for(int v = 0; v < N; ++v)
        z[v] = PoissonBracket(x, y[v], trunc);
    return z;
}


template <class T, int N>
std::istream &operator>>(std::istream &is, FVps<T, N> &vps) {
    return vps.get(is);
}


template <class T, int N>
std::ostream &operator<<(std::ostream &os, const FVps<T, N> &vps) {
    return vps.put(os);
}

#endif // CLASSIC_FVps_CC
