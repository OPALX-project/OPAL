#ifndef CLASSIC_FTps_CC
#define CLASSIC_FTps_CC

// ------------------------------------------------------------------------
// $RCSfile: FTps.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.5.2.15 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template class: FTps<T,N>
//   Representation of a truncated power series
//   with values of type T and a fixed number of variables N.
//
// ------------------------------------------------------------------------
// Class category: FixedAlgebra
// ------------------------------------------------------------------------
//
// $Date: 2003/11/07 18:06:20 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

//#define DEBUG_FTps_CC

#include "FixedAlgebra/FTps.h"
#include "FixedAlgebra/FTpsData.h"
#include "Algebra/Array1D.h"
#include "FixedAlgebra/FArray1D.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FMonomial.h"
#include "FixedAlgebra/FVector.h"
#include "Utilities/ConvergenceError.h"
#include "Utilities/DivideError.h"
#include "Utilities/FormatError.h"
#include "Utilities/LogicalError.h"
#include <climits>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <functional>
#include <string>


// Representation class FTpsRep<T,N>.
// The representation of an FTps<T,N> is based on a mechanism which avoids
// frequent memory allocations.
// This code contains various tweaks for speed.
// ------------------------------------------------------------------------

template<class T, int N>
const int FTps<T, N> :: EXACT = INT_MAX;

template <class T, int N>
class FTpsRep {

    // The only classes allowed access to the representation.
    friend class FTps<T, N>;
    friend class FVps<T, N>;

    // Constructor.
    FTpsRep(int minOrder, int maxOrder, int trcOrder);

    // Destructor.
    ~FTpsRep();

    // Clear FTpsRep<T,N> to zero.
    void clear();

    // Clear FTpsRep<T,N> orders minOrder through maxOrder to zero.
    void clear(int minOrder, int maxOrder);

    // Return beginning of monomial array.
    T *begin() const
    { return data; }

    // Return end of monomial array.
    T *end() const
    { return data + len; }

    // Return beginning of coefficient storage for given order.
    T *begin(int order) const
    { return data + FTpsData<N>::orderStart(order); }

    // Return end of coefficient storage for given order.
    T *end(int order) const
    { return data + FTpsData<N>::orderEnd(order); }

    // Pointer to next object in memory pool, when free.
    FTpsRep<T, N> *next;

    // Reference count.
    int ref;

    // Orders of this object.
    int minOrd;    // The lowest order which contains non-zero coefficients.
    int maxOrd;    // The highest order which contains non-zero coefficients.
    int trcOrd;    // The truncation (precision) order, beyond which terms are truncated.
    int allocOrd;  // The allocation order: the maximum order for which this FTpsRep has memory.
    // Note: We must have: 0 <= minOrd <= maxOrd <= allocOrd < EXACT,
    //                                    maxOrd <= trcOrd,
    //                     and   trcOrd _either_ <= allocOrd _or_ = EXACT.
    //       Also, if trcOrd == EXACT, then (initially) allocOrd == maxOrd.

    // Length of coefficient array.
    int len;

    // Coefficient array.
    T *data;
};


// Implementation of representation class FTpsRep<T,N>.
// ------------------------------------------------------------------------

template <class T, int N>
FTpsRep<T, N>::FTpsRep(int minOrder, int maxOrder, int trcOrder):
    next(0), ref(1), minOrd(minOrder), maxOrd(maxOrder), trcOrd(trcOrder) {
    allocOrd = trcOrd;
    if(trcOrd == FTps<T, N>::EXACT) allocOrd = maxOrd;
    FTpsData<N>::setup(allocOrd);
    len = FTpsData<N>::getSize(allocOrd);
    data = new T[len];
}


template <class T, int N>
FTpsRep<T, N>::~FTpsRep() {
    delete [] data;
}


template <class T, int N>
void FTpsRep<T, N>::clear() {
    clear(minOrd, maxOrd);
}


template <class T, int N>
void FTpsRep<T, N>::clear(int minOrder, int maxOrder) {
    std::fill(data + FTpsData<N>::orderStart(minOrder),
              data + FTpsData<N>::orderEnd(maxOrder), T(0));
}



// Implementation of template class FTps<T,N>.
// ------------------------------------------------------------------------

template <class T, int N>
FTpsRep<T, N> *FTps<T, N>::freeList[100];

template <class T, int N>
int FTps<T, N>::globalTruncOrder;


template <class T, int N>
FTps<T, N>::FTps() {
    itsRep = allocate(0, 0, EXACT);
    itsRep->data[0] = T(0);
}


template <class T, int N>
FTps<T, N>::FTps(const FTps<T, N> &rhs) {
    itsRep = rhs.itsRep;
    itsRep->ref++;
}


template <class T, int N>
FTps<T, N>::FTps(int minOrder, int maxOrder, int trcOrder) {
    itsRep = allocate(minOrder, maxOrder, trcOrder);
    itsRep->clear(minOrder, maxOrder);
}


template <class T, int N>
FTps<T, N>::FTps(const T &rhs) {
    itsRep = allocate(0, 0, EXACT);
    itsRep->data[0] = rhs;
}


template <class T, int N>
FTps<T, N>::FTps(int rhs) {
    itsRep = allocate(0, 0, EXACT);
    itsRep->data[0] = T(rhs);
}


template <class T, int N>
FTps<T, N>::~FTps() {
    itsRep->ref--;
    if(itsRep->ref == 0) deallocate(itsRep);
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator=(const FTps<T, N> &rhs) {
    if(itsRep != rhs.itsRep) {
        itsRep->ref--;
        if(itsRep->ref == 0) deallocate(itsRep);
        itsRep = rhs.itsRep;
        itsRep->ref++;
    }
    return *this;
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator=(const T &rhs) {
    itsRep->ref--;
    if(itsRep->ref == 0) deallocate(itsRep);
    itsRep = allocate(0, 0, EXACT);
    itsRep->data[0] = T(rhs);
    return *this;
}


template <class T, int N>
const T FTps<T, N>::getCoefficient(int index) const {
    if(index < 0) {
        std::cerr << " <*** WARNING ***> from FTps<T,N>::getCoefficient(index):\n"
                  << "    No coefficient has a negative index; returning 0." << std::endl;
        return T(0);
    }

    if(index < orderStart(itsRep->minOrd) || orderEnd(itsRep->maxOrd) <= index) return T(0);
    return itsRep->data[index];
}


template <class T, int N>
void FTps<T, N>::setCoefficient(int index, const T &value) {
    if(index < 0) {
        std::cerr << " <*** WARNING ***> from FTps<T,N>::setCoefficient(index, value):\n"
                  << "    Ignoring request because index < 0." << std::endl;
        return;
    }

    int order = FTpsData<N>::getOrder(index);

    // Ignore requests beyond truncation order.
    if(order > itsRep->trcOrd) return;

    // If necessary, allocate more space.
    // (This happens _only_ if itsRep->trcOrd == EXACT.)
    if(order > itsRep->allocOrd) grow(order, EXACT);

    // Set coefficient.
    unique();
    if(order < itsRep->minOrd) setMinOrder(order);
    else if(order > itsRep->maxOrd) setMaxOrder(order);
    itsRep->data[index] = value;
}


template <class T, int N>
const T FTps<T, N>::getCoefficient(const FMonomial<N> &mono) const {
    int index = FTpsData<N>::getIndex(mono);
    return getCoefficient(index);
}


template <class T, int N>
void FTps<T, N>::setCoefficient(const FMonomial<N> &mono, const T &value) {
    int index = FTpsData<N>::getIndex(mono);
    setCoefficient(index, value);
}


template <class T, int N>
inline const T FTps<T, N>::operator[](int index) const {
    return itsRep->data[index];
}


template <class T, int N>
inline T &FTps<T, N>::operator[](int index) {
    unique();
    return itsRep->data[index];
}


template <class T, int N>
const T FTps<T, N>::operator[](const FMonomial<N> &mono) const {
    int index = FTpsData<N>::getIndex(mono);
    return itsRep->data[index];
}


template <class T, int N>
T &FTps<T, N>::operator[](const FMonomial<N> &mono) {
    unique();
    int index = FTpsData<N>::getIndex(mono);
    return itsRep->data[index];
}


template <class T, int N>
int FTps<T, N>::getSize(int order) {
    return FTpsData<N>::getSize(order);
}


template <class T, int N>
const FMonomial<N> &FTps<T, N>::getExponents(int index) {
    return FTpsData<N>::getExponents(index);
}


template <class T, int N>
int FTps<T, N>::getIndex(const FMonomial<N> &mono) {
    return FTpsData<N>::getIndex(mono);
}


template <class T, int N>
void FTps<T, N>::setMinOrder(int order) {
    // Check sanity.
    if(order == EXACT)
        throw LogicalError("FTps<T,N>::setMinOrder(order)", "Cannot set minimum order to EXACT.");
    if(order < 0)
        throw LogicalError("FTps<T,N>::setMinOrder(order)", "Cannot set a negative minimum order.");
    if(order > globalTruncOrder)
        throw LogicalError("FTps<T,N>::setMinOrder(order)",
                           "Cannot set minimum order beyond globalTruncOrder.");

    unique();
    // Ignore requests beyond truncation order, but issue a warning.
    if(order > itsRep->trcOrd) {
        std::cerr << " <*** WARNING ***> from FTps<T,N>::setMinOrder(order):\n"
                  << "    Cannot set minimum order above truncation order;\n"
                  << "    raising both minimum and maximum orders to truncation order." << std::endl;
        order = itsRep->trcOrd;
    }
    // If necessary, allocate more space.
    // (This can happen _only_ if itsRep->trcOrd == EXACT.)
    else if(order > itsRep->allocOrd) grow(order, EXACT);

    // If raising minOrd above maxOrd, drag max with it and clear that order.
    if(order > itsRep->maxOrd) {
        itsRep->maxOrd = order;
        itsRep->clear(order, order);
    }
    // If decreasing minOrd, stuff zeroes into new coefficients.
    else if(order < itsRep->minOrd) itsRep->clear(order, itsRep->minOrd - 1);

    // Set new minimum order.
    itsRep->minOrd = order;
}


template <class T, int N>
void FTps<T, N>::setMaxOrder(int order) {
    // Check sanity.
    if(order == EXACT)
        throw LogicalError("FTps<T,N>::setMaxOrder(order)", "Cannot set maximum order to EXACT.");
    if(order < 0)
        throw LogicalError("FTps<T,N>::setMaxOrder(order)", "Cannot set a negative maximum order.");
    if(order > globalTruncOrder)
        throw LogicalError("FTps<T,N>::setMaxOrder(order)",
                           "Cannot set maximum order beyond globalTruncOrder.");

    unique();
    // Ignore requests beyond truncation order, but issue a warning.
    if(order > itsRep->trcOrd) {
        std::cerr << " <*** WARNING ***> from FTps<T,N>::setMaxOrder(order):\n"
                  << "    Cannot set maximum order above truncation order;\n"
                  << "    raising maximum order to truncation order." << std::endl;
        order = itsRep->trcOrd;
    }
    // If necessary, allocate more space.
    // (This can happen _only_ if itsRep->trcOrd == EXACT.)
    else if(order > itsRep->allocOrd) grow(order, EXACT);

    // If increasing maxOrd, stuff zeroes into new coefficients.
    if(order > itsRep->maxOrd) itsRep->clear(itsRep->maxOrd + 1, order);
    // If lowering maxOrd below minOrd, drag min with it and clear that order.
    else if(order < itsRep->minOrd) {
        itsRep->minOrd = order;
        itsRep->clear(order, order);
    }

    // Set new maximum order.
    itsRep->maxOrd = order;
}


template <class T, int N>
void FTps<T, N>::setTruncOrder(int order) {
    // Check sanity.
    if(order != EXACT && order < 0)
        throw LogicalError("FTps<T,N>::setTruncOrder(order)", "Cannot set a negative truncation order.");
    if(order != EXACT && order > globalTruncOrder)
        throw LogicalError("FTps<T,N>::setTruncOrder(order)",
                           "Cannot set truncation order beyond globalTruncOrder.");

    unique();
    // If lowering trcOrd below minOrd, drag min and max with it, and clear that order.
    if(order < itsRep->minOrd) {
        itsRep->minOrd = order;
        itsRep->maxOrd = order;
        itsRep->clear(order, order);
    }
    // If lowering trcOrd below maxOrd, drag max with it.
    else if(order < itsRep->maxOrd) itsRep->maxOrd = order;
    // If necessary, allocate more space.
    else if(order != EXACT && order > itsRep->allocOrd) grow(itsRep->maxOrd, order);

    // Set new truncation order.
    itsRep->trcOrd = order;
}


template <class T, int N>
void FTps<T, N>::setGlobalTruncOrder(int order) {
    // Check sanity.
    if(order == EXACT)
        throw LogicalError("FTps<T,N>::setGlobalTruncOrder(order)", "Cannot make globalTruncOrder EXACT.");
    if(order <= 0)
        throw LogicalError("FTps<T,N>::setGlobalTruncOrder(order)",
                           "Cannot set a negative global truncation order.");

    globalTruncOrder = order;
    FTpsData<N>::setup(order);
}


template <class T, int N>
const Array1D<int> &FTps<T, N>::getProductArray(int index) {
    return FTpsData<N>::getProductArray(index);
}


template <class T, int N>
const Array1D<int> &FTps<T, N>::getVariableList(int index) {
    return FTpsData<N>::getVariableList(index);
}


template <class T, int N>
const Array1D<TpsSubstitution> &FTps<T, N>::getSubTable() {
    return FTpsData<N>::getSubTable();
}


template <class T, int N>
FTps<T, N> FTps<T, N>::filter(int minOrder, int maxOrder, int trcOrder) const {
    // Default: trcOrder = EXACT
    checkOrders("filter(minOrder, maxOrder, trcOrder)", minOrder, maxOrder, trcOrder);

    // Compute order limits.
    int myMin = std::max(minOrder, itsRep->minOrd);
    int myMax = std::min(maxOrder, itsRep->maxOrd);
    int myTrc = std::min(trcOrder, itsRep->trcOrd);

    bool OLflag = true;
    if(myMin > myMax) {  // no overlap exists
        OLflag = false;
        myMin = std::min(minOrder, myTrc);
        myMax = std::min(maxOrder, myTrc);
    }
    // Construct filtered FTpsRep.
    FTps<T, N> result(myMin, myMax, myTrc);
    if(OLflag) std::copy(begin(myMin), end(myMax), result.begin(myMin));

    return result;
}


template <class T, int N> inline
FTps<T, N> FTps<T, N>::truncate(int trunc) {
    return filter(0, trunc, trunc);
}


template <class T, int N>
FTps<T, N> FTps<T, N>::makeVariable(int var) {
    FTps<T, N> result(1, 1, EXACT);
    result[var + 1] = T(1);
    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::makeVarPower(int var, int order) {
    FMonomial<N> mono;
    mono[var] = order;
    FTps<T, N> result(order, order, EXACT);
    result[mono] = T(1);
    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::makeMonomial(int index, const T &val) {
    int order = FTpsData<N>::getOrder(index);
    FTps<T, N> result(order, order, EXACT);
    result[index] = val;
    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::makeMonomial(const FMonomial<N> &mono, const T &val) {
    int order = mono.getOrder();
    FTps<T, N> result(order, order, EXACT);
    result[mono] = val;
    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::operator+() const {
    return *this;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::operator-() const {
    int minOrder = getMinOrder();
    int maxOrder = getMaxOrder();
    FTps<T, N> result(minOrder, maxOrder, getTruncOrder());
    std::transform(begin(minOrder), end(maxOrder), result.begin(minOrder), std::negate<T>());
    return result;
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator+=(const FTps<T, N> &rhs) {
    return *this = *this + rhs;
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator-=(const FTps<T, N> &rhs) {
    return *this = *this - rhs;
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator*=(const FTps<T, N> &rhs) {
    return *this = multiply(rhs);
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator/=(const FTps<T, N> &rhs) {
    return *this = divide(rhs);
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator+=(const T &rhs) {
    unique();
    setMinOrder(0);
    itsRep->data[0] += rhs;
    return *this;
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator-=(const T &rhs) {
    unique();
    setMinOrder(0);
    itsRep->data[0] -= rhs;
    return *this;
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator*=(const T &rhs) {
    unique();
    std::transform(begin(getMinOrder()), end(getMaxOrder()), begin(getMinOrder()),
                   std::bind2nd(std::multiplies<T>(), rhs));
    return *this;
}


template <class T, int N>
FTps<T, N> &FTps<T, N>::operator/=(const T &rhs) {
    if(rhs == T(0)) throw DivideError("FTps::operator/=()");
    return *this *= T(1) / rhs;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::scaleMonomials(const FTps<T, N> &rhs) const {
    // Determine orders of result.
    int f_min = getMinOrder();
    int g_min = rhs.getMinOrder();
    int f_max = getMaxOrder();
    int g_max = rhs.getMaxOrder();
    int trcOrder = std::min(getTruncOrder(), rhs.getTruncOrder());
    int minOrder = std::max(f_min, g_min);
    int maxOrder = std::min(f_max, g_max);

    if(minOrder <= maxOrder) {  // overlap exists
        FTps<T, N> result(minOrder, maxOrder, trcOrder);
        std::transform(begin(minOrder), end(maxOrder), rhs.begin(minOrder),
                       result.begin(minOrder), std::multiplies<T>());
        return result;
    } else // no overlap ==> complain big time!
        throw LogicalError("FTps<T,N>::scaleMonomials(rhs)", "No overlap exists.");
}


template <class T, int N>
FTps<T, N> FTps<T, N>::multiplyVariable(int var, int trunc) const {
    // Default: trunc = EXACT

    int f_min = getMinOrder();
    if(trunc <= f_min) {
#ifdef DEBUG_FTps_CC
        std::cerr << " <*** WARNING ***> from FTps<T,N>::multiplyVariable(var, trunc):\n"
                  << "    Result truncated out of existence; returning a zero polynomial."
                  << std::endl;
#endif
        return FTps<T, N>(trunc, trunc, trunc);
    } else {
        // Determine orders of result.
        int trcOrder = getTruncOrder();
        if (trcOrder == EXACT) {
            trcOrder = trunc;
        } else {
            ++ trcOrder;
            trcOrder = std::min(trcOrder, trunc);
        }
        // trcOrder = (trcOrder == EXACT) ? trunc : std::min((++trcOrder), trunc);
        int maxOrder = std::min(1 + getMaxOrder(), trcOrder);

        // Allocate result.
        FTps<T, N> result(1 + f_min, maxOrder, trcOrder);
        const T *f = begin();
        T *g = result.begin();
        const Array1D<int> &product = FTpsData<N>::getProductArray(var + 1);

        // Do multiplication.
        int ie = orderStart(maxOrder);
        for(int i = orderStart(f_min); i < ie; i++)
            g[product[i]] = f[i];
        return result;
    }

}


template <class T, int N>
FTps<T, N> FTps<T, N>::multiply(const FTps<T, N> &rhs, int trunc) const {
    // Default: trunc = EXACT

    //Determine orders of result.
    int trcOrder;
    if(itsRep->trcOrd == EXACT) {
        if(rhs.itsRep->trcOrd == EXACT) trcOrder = EXACT;
        else trcOrder = itsRep->minOrd + rhs.itsRep->trcOrd;
    } else if(rhs.itsRep->trcOrd == EXACT) trcOrder = itsRep->trcOrd + rhs.itsRep->minOrd;
    else trcOrder = std::min(itsRep->trcOrd + rhs.itsRep->minOrd, itsRep->minOrd + rhs.itsRep->trcOrd);
    trcOrder = std::min(trcOrder, trunc);
    int maxOrder = std::min(itsRep->maxOrd + rhs.itsRep->maxOrd, trcOrder);
    int minOrder = itsRep->minOrd + rhs.itsRep->minOrd;

    if(minOrder <= maxOrder) {
        // Allocate result.
        FTps<T, N> result(minOrder, maxOrder, trcOrder);

        // Assign data pointers and starting indices.
        T *f = begin();
        T *g = rhs.begin();
        T *h = result.begin();
        int first_f = orderStart(itsRep->minOrd);
        int first_g = orderStart(rhs.itsRep->minOrd);

        int f_max = itsRep->maxOrd;
        int g_max = std::min(rhs.itsRep->maxOrd, maxOrder);

        // Loop over orders in rhs.
        for(int gOrd = rhs.itsRep->minOrd; gOrd <= g_max; gOrd++) {
            // Set ending indices.
            int last_f = orderEnd(std::min(f_max, maxOrder - gOrd));
            int last_g = orderEnd(gOrd);
            // Do multiplies for each entry (of order gOrd) in rhs.
            for(int j = first_g; j < last_g; j++) {
                T gj = g[j];
                if(gj != T(0)) {
                    const Array1D<int> &prod = FTpsData<N>::getProductArray(j);
                    for(int i = first_f; i < last_f; i++) h[prod[i]] += f[i] * gj;
                }
            }
            // Reset starting index for next order of rhs.
            first_g = last_g;
        }
        return result;
    } else {
#ifdef DEBUG_FTps_CC
        std::cerr << " <*** WARNING ***> from FTps<T,N>::multiply(rhs, trunc):\n"
                  << "    Result truncated out of existence; returning a zero polynomial."
                  << std::endl;
#endif
        return FTps<T, N>(trcOrder, trcOrder, trcOrder);
    }
}


template <class T, int N>
FTps<T, N> FTps<T, N>::inverse(int trunc) const {
    // Default: trunc = EXACT

    // Check sanity.
    T b0 = itsRep->data[0];
    if(b0 == T(0) || itsRep->minOrd != 0) {
        std::cerr << " <*** ERROR ***> in FTps::inverse():\n"
                  << "    Cannot invert a polynomial with zero constant term." << std::endl;
        throw DivideError("FTps::inverse()");
    }

    // Get orders.
    int trcOrder = std::min(itsRep->trcOrd, trunc);
    int maxOrder = std::min(itsRep->maxOrd, trcOrder);
    if(trcOrder > globalTruncOrder)
        throw LogicalError("FTps::inverse(trunc)", "Truncation order exceeds globalTruncOrder!");

    // Allocate result.
    FTps<T, N> result(0, trcOrder, trcOrder);

    // Initialisations.
    const T *b = itsRep->data;
    T *c = result.itsRep->data;
    c[0] = T(1) / b0;
    T neg_c0 = -c[0];

    // Return trivial result.
    if(trcOrder == 0) return result;

    // Loop over orders to construct.
    for(int m = 1; m <= trcOrder; ++m) {
        // Use the c_k's from previous orders to build this order.
        for(int mc = 0; mc < m; ++mc) {
            int mb = m - mc;
            if(mb > maxOrder) continue;
            int start_c = orderStart(mc), end_c = orderEnd(mc);
            int start_b = orderStart(mb), end_b = orderEnd(mb);
            for(int k = start_c; k < end_c; ++k) {
                T ck = c[k];
                if(ck != T(0)) {
                    const Array1D<int> &prod = FTpsData<N>::getProductArray(k);
                    for(int i = end_b; i-- > start_b;) c[prod[i]] += b[i] * ck;
                }
            }
        }
        int je = orderEnd(m);
        for(int j = orderStart(m); j < je; ++j) c[j] *= neg_c0;
    }

    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::divide(const FTps<T, N> &rhs, int trunc) const {
    // Default: trunc = EXACT

    // Algorithm: To compute C = A / B, we may form A = B * C and equate coefficients
    // of like terms on each side.  The resulting system of equations has a "triangular"
    // nature and can be solved explicitly in an order-by-order fashion.  The
    // coefficients of C are given by
    //               c_0 = a_0 / b_0,
    //           and c_j = (a_j - \sum_{i+k=j, i!=0} b_i * c_k) / b_0, for j!=0.
    // These formulae hold for multivariate polynomials if we interpret the subscripts
    // j, i, and k as exponent vectors.

    // Check sanity.
    T b0 = rhs.itsRep->data[0];
    if(b0 == T(0) || rhs.itsRep->minOrd != 0) {
        std::cerr << " <*** ERROR ***> in FTps::divide(rhs,trunc):\n"
                  << "    Cannot divide by a polynomial with zero constant term." << std::endl;
        throw DivideError("FTps::divide()");
    }

    // Get orders.
    int a_min = getMinOrder(), a_max = getMaxOrder(), a_trc = getTruncOrder();
    int b_max = rhs.getMaxOrder(), b_trc = rhs.getTruncOrder();
    int c_trc;
    if(b_trc == EXACT) c_trc = std::min(a_trc, trunc);
    else { c_trc = std::min(a_trc, b_trc + a_min);  c_trc = std::min(c_trc, trunc); }
    if(c_trc > globalTruncOrder)
        throw LogicalError("FTps::divide(rhs,trunc)", "Truncation order exceeds globalTruncOrder!");
    int c_max = c_trc, c_min = a_min;
    int ac_max_min = std::min(a_max, c_max);

    // Allocate result.
    FTps<T, N> result(c_min, c_max, c_trc);

    // Initialisations.
    std::copy(begin(c_min), end(ac_max_min), result.begin(c_min)); // start w/ C = A
    const T *b = rhs.itsRep->data;
    T *c = result.itsRep->data;
    T b0inv = T(1) / b0;
    if(c_min == 0) c[0] *= b0inv;

    // Return trivial result.
    if(c_trc == 0) return result;

    // Loop over orders to construct.
    int mi = std::max(1, c_min);
    for(int m = mi; m <= c_max; ++m) {
        // Use the lower-order c_k's to build the present order.
        for(int mc = c_min; mc < m; ++mc) {
            int mb = m - mc;
            if(mb > b_max) continue;
            int start_c = orderStart(mc), end_c = orderEnd(mc);
            int start_b = orderStart(mb), end_b = orderEnd(mb);
            for(int k = start_c; k < end_c; ++k) {
                T ck = c[k];
                if(ck != T(0)) {
                    const Array1D<int> &prod = FTpsData<N>::getProductArray(k);
                    for(int i = end_b; i-- > start_b;) c[prod[i]] -= b[i] * ck;
                }
            }
        }
        int je = orderEnd(m);
        for(int j = orderStart(m); j < je; ++j) c[j] *= b0inv;
    }

    return result;
}


template <class T, int N>
bool FTps<T, N>::operator==(const FTps<T, N> &rhs) const {
    if(itsRep == rhs.itsRep) return true;
    // else ...

    // Determine ranges (of indices) to test.
    int r_trc = std::min(getTruncOrder(), rhs.getTruncOrder());
    int f_min = orderStart(getMinOrder()),     f_max = orderEnd(std::min(getMaxOrder(), r_trc));
    int g_min = orderStart(rhs.getMinOrder()), g_max = orderEnd(std::min(rhs.getMaxOrder(), r_trc));
    int r_min = std::min(f_min, g_min),         r_max = std::max(f_max, g_max);

    const T *f = begin();
    const T *g = rhs.begin();

    // Test *this == rhs (f == g).
    if(g_min < f_max && f_min < g_max) {  // overlap exists
        int i = r_min, fg_max_min = std::min(f_max, g_max);
        for(; i < g_min; i++)      if(f[i] != T(0)) return false;     //: test below
        for(; i < f_min; i++)      if(g[i] != T(0)) return false;     //     overlap
        for(; i < fg_max_min; i++) if(f[i] != g[i]) return false;     //: test overlap
        for(; i < f_max; i++)      if(f[i] != T(0)) return false;     //: test above
        for(; i < g_max; i++)      if(g[i] != T(0)) return false;     //     overlap
    } else { // gap exists
        if(f_max <= g_min) {  // |f| < |g|
            for(int i = r_min ; i < f_max; i++) if(f[i] != T(0)) return false;
            for(int i = g_min ; i < r_max; i++) if(g[i] != T(0)) return false;
        } else { // |g| < |f|
            for(int i = r_min; i < g_max; i++) if(g[i] != T(0)) return false;
            for(int i = f_min; i < r_max; i++) if(f[i] != T(0)) return false;
        }
    }

    return true;
}


template <class T, int N>
bool FTps<T, N>::operator==(const T &rhs) const {
    int f_min = getMinOrder();
    const T *f = begin();

    if(f_min == 0) {
        if(*f != rhs) return false;
        f_min++;
    } else if(rhs != T(0)) return false;

    int i = orderStart(f_min);
    int iend = orderEnd(getMaxOrder());
    for(; i < iend; i++)
        if(f[i] != T(0)) return false;

    return true;
}


template <class T, int N>
bool FTps<T, N>::operator!=(const FTps<T, N> &rhs) const {
    return !(*this == rhs);
}


template <class T, int N>
bool FTps<T, N>::operator!=(const T &rhs) const {
    return !(*this == rhs);
}


template <class T, int N>
Array1D<T> FTps<T, N>::evalMonoms(const FVector<T, N> &rhs, int maxOrder) {
    int size = getSize(maxOrder);

    // Allocate result.
    Array1D<T> result(size);

    // Initialize pointers to monomials and argument vector.
    T *mon = result.begin();
    const T *z = rhs.begin();

    // Constant term.
    *mon++ = T(1);

    // Linear terms.
    if(maxOrder > 0) std::copy(z, z + N, mon);

    // Higher-order terms.
    if(maxOrder > 1) {
        mon += N;
        int N1 = N - 1;
        int m = 2;
        for(; m <= maxOrder; ++m) {
            // Build monomials of order m from those of order m-1.
            T *monx = mon;                                    // Remember where to stop.
            for(int k = 0, k1 = N1; k < N; k++, k1--) {
                T zk = z[k];                                    // The zk times monomials
                T *mona = monx - FTpsData<N>::orderStart(m, k1); // [mona..monx) yield the
                while(mona != monx) *mon++ = zk * *mona++;      // next ones of order m.
            }
        }
    }

    return result;
}

template <class T, int N>
T FTps<T, N>::evaluate(const FVector<T, N> &rhs) const {
    // Note: We evaluate the linear part separately
    // not out of necessity, but to speed execution.

    // Get coefficient array, f, argument vector, z, and orders.
    const T *f = begin();
    const T *z = rhs.begin();
    int minOrder = getMinOrder();
    int maxOrder = getMaxOrder();

    // Allocate array for monomials and pointer to run along it.
    // Using static local variable to avoid fragmenting memory.
    // Last allocated monoms array is never deleted, but OS
    // will take care of this when program exits.
    static T *monoms = 0;
    static int myOrd = -1;
    if(maxOrder > myOrd) {
        if(monoms) delete [] monoms;
        monoms = new T[getSize(maxOrder)];
        myOrd = maxOrder;
    }
    T *mon = monoms;

    // Initialize monomials and result.
    *mon++ = T(1);
    T result = T(0);
    if(minOrder == 0) {
        result = *f++;
        if(maxOrder == 0) return result;
    } else f++;

    // Load order one monomials.
    std::copy(z, z + N, mon);

    // Accumulate linear contribution.
    T *monx = mon + N;
    if(minOrder <= 1) while(mon != monx) result += *f++ * *mon++;
    else { f += N; mon = monx; }

    // Evaluate higher-order contributions.
    // First build monomials below minOrder.
    int N1 = N - 1;
    int m = 2;
    for(; m < minOrder; m++) {
        // Build monomials of order m from those of order m-1.
        monx = mon;                                       // Remember where to stop.
        for(int k = 0, k1 = N1; k < N; k++, k1--) {
            T zk = z[k];                                    // The zk times monomials
            T *mona = monx - FTpsData<N>::orderStart(m, k1); // [mona..monx)
            while(mona != monx) {                           // yield the next ones
                *mon++ = zk * *mona++;                        // of order m.
                f++;                                          // Skip coefficient.
            }
        }
    }
    // Now add in the higher-order terms.
    for(; m <= maxOrder; m++) {
        // Build monomials of order m from those of order m-1.
        monx = mon;                                       // Remember where to stop.
        for(int k = 0, k1 = N1; k < N; k++, k1--) {
            T zk = z[k];                                    // The zk times monomials
            T *mona = monx - FTpsData<N>::orderStart(m, k1); // [mona..monx)
            while(mona != monx) {                           // yield the next ones
                *mon = zk * *mona++;                          // of order m.
                result += *f++ * *mon++;                      // Accumulate result.
            }
        }
    }

    return result;
}


template <class T, int N>
Array1D<int> FTps<T, N>::getSubstOrders(const FVps<T, N> &rhs, int trunc) const {
    // Default: trunc = EXACT

    // Get orders.
    Array1D<int> ordersL(3), ordersR(3);
    ordersL[0] = getMinOrder(),     ordersL[1] = getMaxOrder(),     ordersL[2] = getTruncOrder();
    ordersR[0] = rhs.getMinOrder(), ordersR[1] = rhs.getMaxOrder(), ordersR[2] = rhs.getTruncOrder();

    Array1D<int> result = getSubstOrders(ordersL, ordersR, trunc);
    return result;
}


template <class T, int N>
Array1D<int> FTps<T, N>::getSubstOrders(Array1D<int> &ordersL, Array1D<int> &ordersR, int trunc) {
    // Default: trunc = EXACT

    // Get orders.
    int f_min = ordersL[0], f_max = ordersL[1], f_trc = ordersL[2];
    int r_min = ordersR[0], r_max = ordersR[1], r_trc = ordersR[2];
    int g_min = f_min * r_min, g_max = f_max * r_max, g_trc;
    if(f_trc == EXACT) {
        if(r_trc == EXACT) {g_trc = (g_max <= trunc) ? EXACT : trunc;}
        else {
            if(f_max == 0) g_trc = EXACT;
            else {
                g_trc = r_trc;
                if(f_min != 0) g_trc += g_min - r_min;
                g_trc = std::min(g_trc, trunc);
            }
        }
    } else { // f_trc != EXACT
        if(r_min != 0) {
            g_trc = std::min((f_trc + 1) * r_min - 1, trunc);
            if(r_trc != EXACT && f_max != 0) {
                int t_trc = r_trc;
                if(f_min != 0) t_trc += g_min - r_min;
                g_trc = std::min(g_trc, t_trc);
            }
        } else { // r_min == 0
#ifdef DEBUG_FTps_CC
            std::cerr << " <*** WARNING ***> from FTps<T,N>::getSubstOrders(ordersL,ordersR,trunc):\n"
                      << "    Incomplete computation of feed-down terms.\n" << std::endl;
#endif
            if(f_max == 0) g_trc = 0;
            else {
                g_trc = std::min(f_trc, r_trc);
                g_trc = std::min(g_trc, trunc);
            }
        }
    }

    g_max = std::min(g_max, g_trc);
    // NB: If trunc is sufficiently small, it can happen that this last line sets
    // g_max (== g_trc) below g_min.  This will mean that g = 0 + O(z^{g_trc+1}).
    // Any routine that calls this one MUST check for this possibility.

    Array1D<int> orders(3);
    orders[0] = g_min;
    orders[1] = g_max;
    orders[2] = g_trc;

    return orders;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::substitute(const FMatrix<T, N, N> &mat, int n) const {
    // Check sanity.
    if(n < 0)
        throw LogicalError("FTps<T,N>::substitute(mat,n)",
                           "Transformation order, n, is negative.");
    else if(n > globalTruncOrder)
        throw LogicalError("FTps<T,N>::substitute(mat,n)",
                           "Transformation order, n, exceeds globalTruncOrder.");

    // Get orders.
    int minOrder = getMinOrder(), maxOrder = getMaxOrder(), trcOrder = getTruncOrder();

    //Allocate result.
    FTps<T, N> result(minOrder, maxOrder, trcOrder);
    std::copy(begin(minOrder), end(maxOrder), result.begin(minOrder));

    // Return trivial cases.
    if(n > trcOrder) {
#ifdef DEBUG_FTps_CC
        std::cerr << " <*** WARNING ***> from FTps<T,N>::substitute(mat,n):\n"
                  << "    Transformation order exceeds truncation order;\n"
                  << "    returning polynomial unchanged." << std::endl;
#endif
        return result;
    }
    if(n == 0 || n < minOrder || maxOrder < n) return result;

    // Allocate working array; use static
    // local memory to avoid fragmentation.
    // NB: Instead of static local memory, one could use the data storage of an FTps
    // ---it's the correct size, it's reference counted, and it's retained on a stack.
    static T *t = 0;
    static int max_n = -1;
    if(n > max_n) {
        if(t) delete [] t;
        t = new T[getSize(n)];
        max_n = n;
    }

    // Initialisations.
    const T *fj = begin(n);
    T *g = result.begin();
    T *t1 = t + 1;
    const Array1D<int> *oldvrbl = 0;
    int start_n = orderStart(n), end_n = orderEnd(n);

    // Clear order n.
    std::fill(g + start_n, g + end_n, T(0));

    // Loop over order n monomials.
    for(int j = start_n; j < end_n; ++j, ++fj) {
        // Skip monomials with coefficient zero.
        if(*fj == T(0)) continue;

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
        std::fill(t + orderStart(ord), t + end_n, T(0));
        // Build the remainder.
        while(ord <= n) {
            // Build next order part of transformed monomial by multiplying
            // the part that is one order lower by the transformed version
            // of the next variable in the variable list.
            int ord1 = ord - 1;
            int start_l = orderStart(ord1), end_l = orderEnd(ord1);
            mv = mat[(*vrbl)[ord1]];  // transformed version of next variable
            for(int k = 0; k < N; k++) {
                T mvk = mv[k];
                if(mvk == T(0)) continue;
                const Array1D<int> &prod = FTpsData<N>::getProductArray(k + 1);
                for(int l = start_l; l < end_l; l++) t[prod[l]] += mvk * t[l];
            }
            ++ord;
        }
        //Increment g by f[j] * transformed monomial.
        for(int i = start_n; i < end_n; i++) g[i] += *fj * t[i];
        // Save variable list for comparison with the next one.
        oldvrbl = vrbl;
    }

    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::substitute(const FMatrix<T, N, N> &mat, int nl, int nh) const {
    // Check sanity.
    if(nl > nh)
        throw LogicalError("FTps<T,N>::substitute(mat,nl,nh)",
                           "Inconsistent transformation orders: nl > nh.");
    if(nl < 0)
        throw LogicalError("FTps<T,N>::substitute(mat,nl,nh)",
                           "Transformation order nl is negative.");
    else if(nh > globalTruncOrder)
        throw LogicalError("FTps<T,N>::substitute(mat,nl,nh)",
                           "Transformation order nh exceeds globalTruncOrder.");

    // Get orders.
    int minOrder = getMinOrder(), maxOrder = getMaxOrder(), trcOrder = getTruncOrder();

    //Allocate result.
    FTps<T, N> result(minOrder, maxOrder, trcOrder);
    std::copy(begin(minOrder), end(maxOrder), result.begin(minOrder));

    if(nh > trcOrder) {
#ifdef DEBUG_FTps_CC
        std::cerr << " <*** WARNING ***> from FTps<T,N>::substitute(mat,nl,nh):\n"
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
    std::fill(result.begin(nl), result.end(nh), T(0));

    // Allocate working arrays; use static
    // local memory to avoid fragmentation.
    static T *t = 0;
    static const T **fp;
    static int max_nh = -1;
    if(nh > max_nh) {
        if(t) {
            delete [] t;
            delete [] fp;
        }
        t = new T[getSize(nh)];
        fp = new const T*[nh+1];
        max_nh = nh;
    }
    T *t1 = t + 1;

    // Initialisations.
    // Array element fp[m] points to the next order m monomial to transform.
    for(int m = nl; m <= nh; ++m) fp[m] = begin(m);
    T *g = result.begin();
    const Array1D<int> *oldvrbl = 0;
    int start_nh = orderStart(nh), end_nh = orderEnd(nh), nh1 = nh - 1, nh2 = nh - 2;

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
        while(n1 <= n2 && *fp[n1] == 0) ++n1;
        while(n2 >  n1 && *fp[n2] == 0) --n2;
        // Skip if all monomials have coefficient zero.
        if(n1 > n2) {
            for(int m = ni; m <= nh; ++m) ++fp[m];
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
        std::fill(t + orderStart(ord), t + end_nh, T(0));

        // Build the remainder.
        while(ord <= nh) {
            // Build next order part of transformed monomial by multiplying
            // the part that is one order lower by the transformed version
            // of the next variable in the variable list.
            int ord1 = ord - 1;
            int start_l = orderStart(ord1), end_l = orderEnd(ord1);
            mv = mat[(*vrbl)[ord1]];  // transformed version of next variable
            for(int k = 0; k < N; k++) {
                T mvk = mv[k];
                if(mvk == T(0)) continue;
                const Array1D<int> &prod = FTpsData<N>::getProductArray(k + 1);
                for(int l = start_l; l < end_l; ++l) t[prod[l]] += mvk * t[l];
            }
            ++ord;
        }
        //Increment g by f[j] * transformed monomial.
        for(int m = n1; m <= n2; ++m) {
            T fj = *fp[m];
            int start_m = orderStart(m), end_m = orderEnd(m);
            for(int i = start_m; i < end_m; i++) g[i] += fj * t[i];
        }
        // Increment pointers in fp[].
        for(int m = ni; m <= nh; ++m) ++fp[m];
        // Save variable list for comparison with the next one.
        oldvrbl = vrbl;
    }

    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::substitute(const FMatrix<T, N, N> &mat) const {
    return substitute(mat, getMinOrder(), getMaxOrder());
}


template <class T, int N>
FTps<T, N> FTps<T, N>::substitute(const FVps<T, N> &rhs, int trunc) const {
    // Default: trunc = EXACT

    // Get orders.
    int f_min = getMinOrder(), f_max = getMaxOrder();
    Array1D<int> orders = getSubstOrders(rhs, trunc);
    int g_min = orders[0], g_max = orders[1], g_trc = orders[2];

    // Make sure we don't trip over globalTruncOrder.
    if(g_trc != EXACT && g_trc > globalTruncOrder)
        throw LogicalError("FTps::substitute(FVps rhs, int trunc)",
                           "Truncation order exceeds globalTruncOrder!");

    // Return trivial case.
    if(g_min > g_max) return FTps<T, N>(g_trc, g_trc, g_trc);

    //Allocate result.
    FTps<T, N> result(g_min, g_max, g_trc);
    if(f_min == 0) result[0] = *begin();
    if(f_max == 0) return result;

    // Set actual range of orders to transform.
    int nl = f_min, nh = f_max;
    if(nl == 0) nl = 1;

    // Allocate working arrays.
    const T *fp[nh+1];
    Array1D< FTps<T, N> > t(nh + 1);

    // Initialisations.
    // Array element fp[m] points to the next order m monomial to transform.
    for(int m = nl; m <= nh; ++m) fp[m] = begin(m);
    const Array1D<int> *oldvrbl = 0;
    int start_nh = orderStart(nh), end_nh = orderEnd(nh), nh1 = nh - 1, nh2 = nh - 2;

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
        while(n1 <= n2 && *fp[n1] == 0) ++n1;
        while(n2 >  n1 && *fp[n2] == 0) --n2;
        // Skip if all monomials have coefficient zero.
        if(n1 > n2) {
            for(int m = ni; m <= nh; ++m) ++fp[m];
            continue;
        }

        int ord;
        // If vk = 0, we must start at the beginning; otherwise,
        // we may keep the first vk orders stored in t.
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
        //Increment result by f[j] * transformed monomial.
        for(int m = n1; m <= n2; ++m) result += *fp[m] * t[m];
        // Increment pointers in fp[].
        for(int m = ni; m <= nh; ++m) ++fp[m];
        // Save variable list for comparison with the next one.
        oldvrbl = vrbl;
    }

    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::derivative(int var) const {
    // Get orders.
    int f_min = getMinOrder(), f_max = getMaxOrder(), f_trc = getTruncOrder();
    int df_min = (f_min == 0) ? f_min : f_min - 1;
    int df_max = (f_max == 0) ? f_max : f_max - 1;
    int df_trc = (f_trc == 0 || f_trc == EXACT) ? f_trc : f_trc - 1;
    if(f_min == 0 && f_min < f_max) f_min = 1;

    //Allocate result.
    FTps<T, N> result(df_min, df_max, df_trc);
    // Return trivial case.
    if(f_max == 0) return result;

    // Initialisations.
    const int *product = FTpsData<N>::getProductArray(var + 1).begin();
    const T *f = begin();
    T *df = result.begin();
    T *dfe = df + orderEnd(df_max);
    int ks = orderStart(df_min);
    df += ks;
    product += ks;

    // Compute derivative.
    // For the contribution to the k-th monomial in the derivative, we look
    // at the monomial whose index is the k-th entry in the product array for
    // variable var.  The var-th entry in that monomial's exponent array is
    // the power brought down by the differentiation.
    while(df != dfe) {
        int kp = *product++;
        *df++ = T(FTpsData<N>::getExponents(kp)[var]) * *(f + kp);
    }

    return result;
}


template <class T, int N>
FVps<T, N> FTps<T, N>::gradient() const {
    FVps<T, N> result;
    for(int i = 0; i < N; ++i) result[i] = derivative(i);
    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::integral(int var, int trunc) const {
    // Default: trunc = EXACT

    // Get orders.
    int f_min = getMinOrder(), f_max = getMaxOrder(), f_trc = getTruncOrder();
    int minOrder = std::min(f_min + 1, trunc);
    int maxOrder = std::min(f_max + 1, trunc);
    int trcOrder;
    if(f_trc == EXACT) trcOrder = trunc;
    else trcOrder = std::min(f_trc + 1, trunc);
    if(maxOrder > globalTruncOrder || (trcOrder != EXACT && trcOrder > globalTruncOrder))
        throw LogicalError("FTps::integral(int var, int trunc)", "Some order exceeds globalTruncOrder!");

    // Allocate result.
    FTps<T, N> result(minOrder, maxOrder, trcOrder);

    // Return trivial result.
    if(f_min >= maxOrder) return result;

    // Initialisations.
    const T *f = begin();
    T *r = result.begin();
    int is = orderStart(f_min);
    int ie = orderEnd(maxOrder - 1);

    // Compute integral.
    // The contribution from the k-th monomial we insert into the integral
    // at the monomial whose index is the k-th entry in the product array
    // for variable var.  The var-th entry in that (the result) monomial's
    // exponent array is the power kicked up by the integration.
    const Array1D<int> &product = FTpsData<N>::getProductArray(var + 1);
    for(int i = is; i < ie; ++i) {
        int k = product[i];
        r[k] = f[i] / double(FTpsData<N>::getExponents(k)[var]);
    }
    return result;
}


template <class T, int N>
FTps<T, N> FTps<T, N>::taylor(const Array1D<T> &series, int order) const {
    FTps<T, N> x(*this);
    x.unique();
    x.itsRep->data[0] = T(0);
    FTps<T, N> result(series[order]);

    for(int m = 1; m <= order; ++m) {
        if(result.itsRep->trcOrd < m) result.setTruncOrder(m);
        result = x.multiply(result, m);

        /*
         * has to be set, otherwise if coefficient is zero --> next coefficient doesn't get
         * counted, e.g. with cos(x)
         */
        result.setMinOrder(0);

        result.itsRep->data[0] = series[order-m];
    }

    return result;
}


template <class T, int N> inline
void FTps<T, N>::unique() {
    if(itsRep->ref > 1) {
        FTpsRep<T, N> *oldRep = itsRep;
        oldRep->ref--;
        itsRep = FTps<T, N>::allocate(itsRep->minOrd, itsRep->maxOrd, itsRep->trcOrd);
        std::copy(oldRep->begin(), oldRep->end(itsRep->maxOrd), itsRep->begin());
    }
}

template <class T, int N>
std::list<int> FTps<T, N>::getListOfNonzeroCoefficients() const {

    // get total number of coefficients
    int size = getSize();

    // initialize list
    std::list<int> coeffs;

    // loop over all coefficients
    for (int i = 0; i < size; ++i) {
        // get index of non-zero coefficients
        if (getCoefficient(i) != 0) {
            coeffs.push_back(i);
        }
    }

    return coeffs;
}

template <class T, int N>
FArray1D<int, N> FTps<T, N>::extractExponents(int index) const {

    // check index
    if ( index < 0 || getSize() - 1 < index)
        throw LogicalError("FVps<T,N>::extractExponents(var)","Index out of range.");

    // get exponents of monomial
    FMonomial<N> mono = FTps<T, N>::getExponents(index);

    // array of exponents
    FArray1D<int, N> expons;

    // copy monomials to array
    for (int i = 0; i < N; ++i)
        expons[i] = mono[i];

    return expons;
}

template <class T, int N>
FTps<T, N> FTps<T, N>::makePower(int power) const {

    if (power < 0)
        throw LogicalError("FTps<T,N>::makePower(power)","Power is negative.");

    FTps<T, N> result = *this;

    // gets truncated in case of being bigger than global truncation order
    for (int i = 1; i < power; ++i) {
        //result *= *this;
        result = result.multiply(*this, globalTruncOrder);
    }
    return result;
}

template <class T, int N>
std::istream &FTps<T, N>::get(std::istream &is) {
    //  is.flags(std::ios::skipws);
    //  char head[4];
    //  is.get(head, 4);
    //
    //  // Check for Tps.
    //  if (strcmp(head, "Tps") != 0)
    //    throw FormatError("FTps::get()", "Flag word \"Tps\" missing.");

    // Check for Tps.
    std::string head;
    is >> head;
    if(head != "Tps")
        throw FormatError("FTps::get()", "Flag word \"Tps\" missing.");

    // Read in order data.
    int minOrder, maxOrder, trcOrder, nVar;
    std::string trunc;
    is >> minOrder >> maxOrder >> trunc >> nVar;
    if(trunc == "EXACT") trcOrder = EXACT;
    else trcOrder = std::atoi(trunc.c_str());

    // Check sanity.
    if(nVar != N)
        throw FormatError("FTps::get()", "Invalid number of variables.");

    // Allocate result.
    FTps<T, N> result(minOrder, maxOrder, trcOrder);
    T coeff;
    FMonomial<N> mono(nVar);
    bool done = false;
    bool fail = false;
    bool first = true;

    while(true) {
        is >> coeff;
        fail = is.fail();

        int order = 0;
        for(int var = 0; var < nVar; var++) {
            int p;
            is >> p;
            fail |= is.fail();
            if(p < 0) done = true;
            mono[var] = p;
            order += mono[var];
        }

        if(fail) throw FormatError("FTps::get()", "File read error");
        if(done) break;

        if(coeff != T(0)) {
            if(first) {
                first = false;
                minOrder = maxOrder = order;
            } else {
                if(order < minOrder) minOrder = order;
                else if(order > maxOrder) maxOrder = order;
            }
            int index = FTpsData<N>::getIndex(mono);
            result[index] = coeff;
        }
    }

    if(minOrder < result.getMinOrder())
        throw FormatError("FTps::get()", "minOrder incorrect in file");
    if(maxOrder > result.getMaxOrder())
        throw FormatError("FTps::get()", "maxOrder incorrect in file");
    result.itsRep->minOrd = minOrder;
    result.itsRep->maxOrd = maxOrder;

    *this = result;
    return is;
}


template <class T, int N>
std::ostream &FTps<T, N>::put(std::ostream &os) const {
    // Print out order data.
    os << "Tps " << itsRep->minOrd << ' ' << itsRep->maxOrd << ' ';
    if(itsRep->trcOrd == EXACT)
        os << "EXACT";
    else
        os << itsRep->trcOrd;
    os << "  " << N << std::endl;
    //  os << "<*** allocOrd = " << itsRep->allocOrd << " ***>" << std::endl;

    // Save old and set new formats.
    std::streamsize old_prec = os.precision(14);
    os.setf(std::ios::scientific, std::ios::floatfield);

    // Print out FTps coefficients and exponents.
    int i = orderStart(itsRep->minOrd);
    const T *f = begin(getMinOrder());
    const T *fe = end(getMaxOrder());
    for(; f != fe; f++, i++) {
        if(*f != T(0)) {
            os << std::setw(24) << *f;
            for(int var = 0; var < N; var++)
                os << std::setw(3) << FTpsData<N>::getExponents(i)[var];
            os << std::endl;
        }
    }
    // End FTps w/ "0.00 -1 .. -1".
    os << std::setw(24) << T(0);
    for(int var = 0; var < N; var++) os << std::setw(3) << (-1);
    os << std::endl;

    // Restore old formats.
    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);

    return os;
}


// Template class FTps<T,N>, private methods.
// ------------------------------------------------------------------------

template <class T, int N> inline
FTpsRep<T, N> *
FTps<T, N>::allocate(int minOrder, int maxOrder, int trcOrder) {
    checkOrders("FTps<T,N>::allocate(minOrder, maxOrder, trcOrder)", minOrder, maxOrder, trcOrder);

    // Determine allocation order.
    int allocOrder = trcOrder == EXACT ? maxOrder : trcOrder;

    // Do the allocating.
    FTpsRep<T, N> *rep = freeList[allocOrder];
    if(rep) {
        // There exists a free object of the correct order: take 'rep' off
        // the stack, point freeList to next element in the linked list, ...
        freeList[allocOrder] = rep->next;
        // ... and set the values of next, ref, and *Ord.
        rep->next = 0;
        rep->ref = 1;
        rep->minOrd = minOrder;
        rep->maxOrd = maxOrder;
        rep->trcOrd = trcOrder;
        //rep->allocOrd = allocOrder;  // This should already be the case! (Don't alter rep->allocOrd!)
#ifdef DEBUG_FTps_CC
        if(allocOrder != rep->allocOrd)
            std::cerr << " <*** allocOrder error ***> " << allocOrder << " != " << rep->allocOrd << std::endl;
#endif
    } else {
        // We must allocate a new object.
        rep = new FTpsRep<T, N>(minOrder, maxOrder, trcOrder);
    }

    // The data array for the returned object is uninitialised.
    return rep;
}


template <class T, int N> inline
void FTps<T, N>::deallocate(FTpsRep<T, N> *rep) {
    int order = rep->allocOrd;
    rep->next = FTps<T, N>::freeList[order];
    FTps<T, N>::freeList[order] = rep;
}


template <class T, int N>
void FTps<T, N>::grow(int maxOrder, int trcOrder) {
    int minOrder = itsRep->minOrd;

    // Allocate new representation.
    FTpsRep<T, N> *newRep = allocate(minOrder, maxOrder, trcOrder);

    // Copy old representation.
    T *f  = itsRep->begin(minOrder);
    T *fe = itsRep->end(std::min(itsRep->maxOrd, maxOrder));
    T *g  = newRep->begin(minOrder);
    T *ge = newRep->end(maxOrder);
    while(f < fe) *g++ = *f++;
    while(g < ge) *g++ = T(0);

    // Discard old representation and attach new one.
    itsRep->ref--;
    if(itsRep->ref == 0) deallocate(itsRep);
    itsRep = newRep;
}


template <class T, int N>
Array1D<int> FTps<T, N>::getRepOrders() const {
    Array1D<int> result(4);
    result[0] = itsRep->minOrd;
    result[1] = itsRep->maxOrd;
    result[2] = itsRep->trcOrd;
    result[3] = itsRep->allocOrd;
    return result;
}


template <class T, int N>
void FTps<T, N>::checkOrders(const std::string &method, int minOrder, int maxOrder, int &trcOrder) {
    bool errF = false;
    std::string message;

    // Check that min-, max-, and trcOrder's have the correct relationships.
    if(!(0 <= minOrder && minOrder <= maxOrder && maxOrder <= trcOrder && maxOrder < EXACT)) {
        errF = true;
        message = "Invalid, or inconsistent, arguments.";
    }
    // Make sure we don't trip over globalTruncOrder.
    else if(maxOrder > globalTruncOrder) {
        errF = true;
        message = "The argument maxOrder exceeds globalTruncOrder.";
    } else if(trcOrder != EXACT && trcOrder > globalTruncOrder) {
        errF = true;
        message = "The argument trcOrder has a finite value that exceeds globalTruncOrder.";
    }
    if(errF) {
        std::cerr << "Method " << method << ", called with arguments (";
        if(minOrder == EXACT) std::cerr << "EXACT,";
        else std::cerr << minOrder << ",";
        if(maxOrder == EXACT) std::cerr << "EXACT,";
        else std::cerr << maxOrder << ",";
        if(trcOrder == EXACT) std::cerr << "EXACT";
        else std::cerr << trcOrder;
        std::cerr << ")." << std::endl;
        throw LogicalError(method, message);
    }
}


// Global functions acting on FTps<T> objects.
// ------------------------------------------------------------------------

template <class T, int N>
FTps<T, N> operator+(const FTps<T, N> &lhs, const FTps<T, N> &rhs) {
    // Determine orders of result.
    int f_min = lhs.getMinOrder(), f_max = lhs.getMaxOrder();
    int g_min = rhs.getMinOrder(), g_max = rhs.getMaxOrder();
    int h_trc = std::min(lhs.getTruncOrder(), rhs.getTruncOrder());
    int h_min = std::min(f_min, g_min), h_max = std::min(std::max(f_max, g_max), h_trc);
    if(h_trc < f_min) h_max = g_max;
    else if(h_trc < g_min) h_max = f_max;

    // Allocate result.
    FTps<T, N> result(h_min, h_max, h_trc);

    // Change from orders to indices.
    f_min = FTps<T, N>::orderStart(f_min);
    f_max = FTps<T, N>::orderEnd(std::min(f_max, h_max));
    g_min = FTps<T, N>::orderStart(g_min);
    g_max = FTps<T, N>::orderEnd(std::min(g_max, h_max));
    h_min = std::min(f_min, g_min);
    h_max = std::max(f_max, g_max);

    const T *f = lhs.begin();
    const T *g = rhs.begin();
    T *h = result.begin();

    // Do the "addition" h = f + g.
    int i = h_min;
    if(g_min < f_max && f_min < g_max) {  // overlap exists
        int fg_max_min = std::min(f_max, g_max);
        for(; i < g_min; i++)      h[i] = f[i];          //: copy below
        for(; i < f_min; i++)      h[i] = g[i];          //     overlap
        for(; i < fg_max_min; i++) h[i] = f[i] + g[i];   //: do f + g
        for(; i < f_max; i++)      h[i] = f[i];          //: copy above
        for(; i < g_max; i++)      h[i] = g[i];          //     overlap
    } else { // gap exists
        if(f_max <= g_min) {  // |f| < |g|
            // Check for special case (corresponds to h_trc in gap).
            if(g_min > h_max) g_min = h_max;
            // Then load h.
            for(; i < f_max; i++) h[i] = f[i];
            for(; i < g_min; i++) h[i] = 0;
            for(; i < h_max; i++) h[i] = g[i];
        } else { // |g| < |f|
            // Check for special case (corresponds to h_trc in gap).
            if(f_min > h_max) f_min = h_max;
            // Then load h.
            for(; i < g_max; i++) h[i] = g[i];
            for(; i < f_min; i++) h[i] = 0;
            for(; i < h_max; i++) h[i] = f[i];
        }
    }

    return result;
}


template <class T, int N>
FTps<T, N> operator-(const FTps<T, N> &lhs, const FTps<T, N> &rhs) {
    // Determine orders of result.
    int f_min = lhs.getMinOrder(), f_max = lhs.getMaxOrder();
    int g_min = rhs.getMinOrder(), g_max = rhs.getMaxOrder();
    int h_trc = std::min(lhs.getTruncOrder(), rhs.getTruncOrder());
    int h_min = std::min(f_min, g_min), h_max = std::min(std::max(f_max, g_max), h_trc);
    if(h_trc < f_min) h_max = g_max;
    else if(h_trc < g_min) h_max = f_max;

    // Allocate result.
    FTps<T, N> result(h_min, h_max, h_trc);

    // Change from orders to indices.
    f_min = FTps<T, N>::orderStart(f_min);
    f_max = FTps<T, N>::orderEnd(std::min(f_max, h_max));
    g_min = FTps<T, N>::orderStart(g_min);
    g_max = FTps<T, N>::orderEnd(std::min(g_max, h_max));
    h_min = std::min(f_min, g_min);
    h_max = std::max(f_max, g_max);

    const T *f = lhs.begin();
    const T *g = rhs.begin();
    T *h = result.begin();

    // Do the "subtraction" h = f - g.
    int i = h_min;
    if(g_min < f_max && f_min < g_max) {  // overlap exists
        int fg_max_min = std::min(f_max, g_max);
        for(; i < g_min; i++)      h[i] =  f[i];          //: copy below
        for(; i < f_min; i++)      h[i] = -g[i];          //     overlap
        for(; i < fg_max_min; i++) h[i] =  f[i] - g[i];   //: do f - g
        for(; i < f_max; i++)      h[i] =  f[i];          //: copy above
        for(; i < g_max; i++)      h[i] = -g[i];          //     overlap
    } else { // gap exists
        if(f_max <= g_min) {  // |f| < |g|
            // Check for special case (corresponds to h_trc in gap).
            if(g_min > h_max) g_min = h_max;
            // Then load h.
            for(; i < f_max; i++) h[i] =  f[i];
            for(; i < g_min; i++) h[i] =  0;
            for(; i < h_max; i++) h[i] = -g[i];
        } else { // |g| < |f|
            // Check for special case (corresponds to h_trc in gap).
            if(f_min > h_max) f_min = h_max;
            // Then load h.
            for(; i < g_max; i++) h[i] = -g[i];
            for(; i < f_min; i++) h[i] =  0;
            for(; i < h_max; i++) h[i] =  f[i];
        }
    }

    return result;
}


template <class T, int N>
FTps<T, N> operator+(const FTps<T, N> &lhs, const T &rhs) {
    FTps<T, N> result(lhs);
    return result += rhs;
}


template <class T, int N>
FTps<T, N> operator-(const FTps<T, N> &lhs, const T &rhs) {
    FTps<T, N> result(lhs);
    return result -= rhs;
}


template <class T, int N>
FTps<T, N> operator+(const T &lhs, const FTps<T, N> &rhs) {
    FTps<T, N> result(rhs);
    return result += lhs;
}


template <class T, int N>
FTps<T, N> operator-(const T &lhs, const FTps<T, N> &rhs) {
    FTps<T, N> result(- rhs);
    return result += lhs;
}


template <class T, int N>
FTps<T, N> operator*(const FTps<T, N> &lhs, const FTps<T, N> &rhs) {
    return lhs.multiply(rhs);
    //  return lhs.multiply(rhs, FTps<T,N>::getGlobalTruncOrder());
}


template <class T, int N>
FTps<T, N> operator/(const FTps<T, N> &lhs, const FTps<T, N> &rhs) {
    return lhs.divide(rhs);
}


template <class T, int N>
FTps<T, N> operator*(const FTps<T, N> &lhs, const T &rhs) {
    FTps<T, N> result(lhs);
    return result *= rhs;
}


template <class T, int N>
FTps<T, N> operator/(const FTps<T, N> &lhs, const T &rhs) {
    FTps<T, N> result(lhs);
    return result /= rhs;
}


template <class T, int N>
FTps<T, N> operator*(const T &lhs, const FTps<T, N> &rhs) {
    FTps<T, N> result(rhs);
    return result *= lhs;
}


template <class T, int N>
FTps<T, N> operator/(const T &lhs, const FTps<T, N> &rhs) {
    return rhs.inverse() * lhs;
}


template <class T, int N>
bool operator==(const T &lhs, const FTps<T, N> &rhs) {
    return rhs == lhs;
}


template <class T, int N>
bool operator!=(const T &lhs, const FTps<T, N> &rhs) {
    return rhs != lhs;
}


template <class T, int N> FVps<T, N>
ExpMap(const FTps<T, N> &H, int trunc) {
    return ExpMap(H, FVps<T, N>(), trunc);
}


template <class T, int N> FTps<T, N>
ExpMap(const FTps<T, N> &H, const FTps<T, N> &f, int trunc) {
    // Default: trunc = EXACT

    // Limit number of iterations.
    const int MAX_ITER = 100;

    // We really ought to throw an exception if H contains linear terms and is not exact,
    // but we're just going to complain!!
    bool FD = false;
    if(H.getTruncOrder() != FTps<T, N>::EXACT && H.getMinOrder() < 2) {
        FD = true;
#ifdef DEBUG_FTps_CC
        std::cerr << " <*** WARNING ***> from ExpMap(H,f,trunc):\n"
                  << "    Incomplete computation of feed-down terms.\n" << std::endl;
#endif
    }
    int fd_trc = std::min(H.getTruncOrder() - 1, f.getTruncOrder());

    // Construct dH = grad(H).J, s.t :H:f = dH.grad(f).
    FVps<T, N> dH;
    for(int i = 0; i < N; i += 2) {
        dH[i]   = - H.derivative(i + 1);
        dH[i+1] =   H.derivative(i);
    }

    // Allocate and initialize result.
    FTps<T, N> expHf = f;

    // Initialize variables.
    FTps<T, N> dHkf  = f;
    FTps<T, N> old   = T(0);

    // Compute series; quit computation if we added nothing last time through loop.
    for(int k = 1; expHf != old; ++k) {
        if(k > MAX_ITER)
            throw ConvergenceError("ExpMap(const FTps<T,N> &H, const FTps<T,N> &f)",
                                   "No convergence in ExpMap(H,f)");

        // Don't initialize dHk1f to 0, as that sets minOrder to 0!
        old = expHf;
        FTps<T, N> ddHkf = dHkf.derivative(0);
        if(FD) ddHkf.setTruncOrder(fd_trc);
        FTps<T, N> dHk1f = dH[0].multiply(ddHkf, trunc);
        for(int v = 1; v < N; ++v) {
            ddHkf = dHkf.derivative(v);
            if(FD) ddHkf.setTruncOrder(fd_trc);
            dHk1f += dH[v].multiply(ddHkf, trunc);
            //dHk1f += dH[v].multiply(dHkf.derivative(v),trunc);
        }
        dHkf = dHk1f / T(k); // :H:^{k}/(k!)
        expHf += dHkf;
    }

    return expHf;
}


template <class T, int N>
FTps<T, N> PoissonBracket(const FTps<T, N> &f, const FTps<T, N> &g, int trunc) {
    // Default: trunc = EXACT

    // Determine orders of result.
    int f_min = f.getMinOrder(), f_max = f.getMaxOrder(), f_trc = f.getTruncOrder();
    int g_min = g.getMinOrder(), g_max = g.getMaxOrder(), g_trc = g.getTruncOrder();
    int h_min = std::max(f_min + g_min - 2, 0), h_max = std::max(f_max + g_max - 2, 0), h_trc;
    if(f_trc == FTps<T, N>::EXACT) {
        if(g_trc == FTps<T, N>::EXACT) {h_trc = (h_max <= trunc) ? FTps<T, N>::EXACT : trunc;}
        else h_trc = std::min(std::max(f_min + g_trc - 2, 0), trunc);
    } else if(g_trc == FTps<T, N>::EXACT) h_trc = std::min(std::max(f_trc + g_min - 2, 0), trunc);
    else {
        h_trc = std::min(f_trc + g_min - 2, f_min + g_trc - 2);
        h_trc = std::min(std::max(h_trc, 0), trunc);
    }
    h_max = std::min(h_max, h_trc);

    // Return trivial result.
    if(h_min > h_max) return FTps<T, N>();

    // Allocate result.
    FTps<T, N> h(h_min, h_max, h_trc);

    for(int q = 0; q < N; q += 2) {
        int p = q + 1;
        //h += f.derivative(q) * g.derivative(p) - f.derivative(p) * g.derivative(q);
        h += f.derivative(q).multiply(g.derivative(p), h_trc)
             - f.derivative(p).multiply(g.derivative(q), h_trc);
    }
    return h;
}


template <class T, int N>
std::istream &operator>>(std::istream &is, FTps<T, N> &tps) {
    return tps.get(is);
}


template <class T, int N>
std::ostream &operator<<(std::ostream &os, const FTps<T, N> &tps) {
    return tps.put(os);
}

#endif // CLASSIC_FTps_CC
