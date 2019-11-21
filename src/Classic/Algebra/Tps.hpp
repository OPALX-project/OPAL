#ifndef CLASSIC_Tps_CC
#define CLASSIC_Tps_CC 1

// ------------------------------------------------------------------------
// $RCSfile: Tps.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Template Class: Tps<class T>
//   Truncated power series in n variables of type T.
//   This file contains the template implementation only.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2001/11/27 23:33:32 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include "Algebra/Tps.h"
#include "Algebra/Array1D.h"
#include "Algebra/Matrix.h"
#include "Algebra/TpsData.h"
#include "Algebra/TpsMonomial.h"
#include "Algebra/TpsSubstitution.h"
#include "Algebra/VpsMap.h"
#include "Utilities/DivideError.h"
#include "Utilities/FormatError.h"
#include "Utilities/LogicalError.h"
#include "Utilities/CLRangeError.h"
#include "Utilities/SizeError.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <new>
#include <cstring>
#include <functional>


// Helper functions.
// ------------------------------------------------------------------------

template <class T> inline
void Tps<T>::setMaxOrder(int order) {
    rep->maxOrd = order;
}


template <class T> inline
void Tps<T>::unique() {
    if(rep->ref > 1) {
        rep->ref--;
        rep = rep->clone();
    }
}


// Template class TpsRep<T>.
//
// The representation of a Tps<T> is based on a mechanism which avoids
// double memory allocation (one for the TpsRep<T> object, one for the
// table of monomials). If this causes problems, it can be reverted to
// using a double allocation by adding a data member "T *dat" to TpsRep<T>,
// changing "TpsRep<T>::data()" to return "dat", and adapting the methods
// "TpsRep<T>::operator new()" and "TpsRep<T>::operator delete()" so as
// to allocate/deallocate an array of "T", pointed to by "dat".
//
// A special memory allocator may help for speed.
// ------------------------------------------------------------------------

template <class T> class TpsRep {

    friend class Tps<T>;

    // Memory management.
    static void *operator new(size_t s, size_t extra);
    static void operator delete(void *);
    static TpsRep<T> *create(int maxOrder, int trcOrder, int variables);
    static TpsRep<T> *zero();
    static void release(TpsRep<T> *);

    // Make a copy of the representation.
    TpsRep<T> *clone();

    // Grab a new reference.
    TpsRep<T> *grab();

    // Return the monomial array.
    T *data();

    // The reference count.
    int ref;

    // Order and size of this object.
    int maxOrd;
    int trcOrd;

    // The length of the monomial table.
    int len;

    // The data structure for bookkeeping.
    TpsData *help;

    // Not implemented.
    TpsRep<T> &operator=(const TpsRep<T> &);

    // May need to add this:
    //   T *dat;
};


template <class T> inline
T *TpsRep<T>::data() {
    return reinterpret_cast<T *>(this + 1);
    // Could be changed to:
    //   return dat;
}


template <class T> inline
void *TpsRep<T>::operator new(size_t s, size_t extra) {
    return new char[s + extra * sizeof(double)];
    // Could be changed to:
    //   TpsRep<T> *p = new char[s];
    //   p->dat = new T[extra];
    //   return p;
}


template <class T> inline
void TpsRep<T>::operator delete(void *p) {
    // May have to add this:
    //   delete [] reinterpret_cast<TpsRep<T>*>(p)->dat;
    delete [] reinterpret_cast<char *>(p);
}


template <class T> inline
TpsRep<T> *TpsRep<T>::create(int maxOrder, int trcOrder, int variables) {
    // Construct descriptor and size.
    TpsData *d = 0;
    int s = 1;
    if(variables) {
        d = TpsData::getTpsData(maxOrder, variables);
        s = d->getSize(maxOrder);
    }

    // Allocate representation and fill in data.
    TpsRep<T> *p = new(s) TpsRep<T>;
    p->ref = 1;
    p->maxOrd = maxOrder;
    p->trcOrd = trcOrder;
    p->len = s;
    p->help = d;

    // Fill monomial coefficients with zeroes.
    std::fill(p->data() + 0, p->data() + p->len, T(0));
    return p;
}


template <class T> inline
TpsRep<T> *TpsRep<T>::zero() {
    // Allocate representation and fill in data.
    TpsRep<T> *p = new(1) TpsRep<T>;
    p->ref = 1;
    p->maxOrd = 0;
    p->trcOrd = Tps<T>::EXACT;
    p->len = 1;
    p->help = 0;

    // Fill monomial coefficients with zeroes.
    new(&p->data()[0]) T(0);
    return p;
}


template <class T> inline
TpsRep<T> *TpsRep<T>::clone() {
    // Allocate copy and copy monomial coefficients.
    TpsRep<T> *p = new(len) TpsRep<T>;
    for(int i = 0; i < len; ++i) {
        new(&p->data()[i]) T(data()[i]);
    }

    // Copy limits and descriptor.
    p->ref = 1;
    p->maxOrd = maxOrd;
    p->trcOrd = trcOrd;
    p->len = len;
    p->help = help;
    return p;
}


template <class T> inline
TpsRep<T> *TpsRep<T>::grab() {
    ++ref;
    return this;
}


template <class T> inline
void TpsRep<T>::release(TpsRep<T> *p) {
    if(--(p->ref) <= 0) delete p;
}


// Template class Tps<T>.
// ------------------------------------------------------------------------

template <class T> int Tps<T>::truncOrder = EXACT;


template <class T>
Tps<T>::Tps():
    rep(TpsRep<T>::zero())
{}


template <class T>
Tps<T>::Tps(int maxOrder, int nVar):
    rep(TpsRep<T>::create(maxOrder, EXACT, nVar))
{}


template <class T>
Tps<T>::Tps(const Tps<T> &rhs):
    rep(rhs.rep->grab())
{}


template <class T>
Tps<T>::Tps(const T &rhs):
    rep(TpsRep<T>::zero()) {
    rep->data()[0] = rhs;
}


template <class T>
Tps<T>::Tps(int rhs):
    rep(TpsRep<T>::zero()) {
    rep->data()[0] = T(rhs);
}


template <class T>
Tps<T>::~Tps() {
    TpsRep<T>::release(rep);
}


template <class T>
Tps<T> &Tps<T>::operator=(const Tps<T> &rhs) {
    if(rep != rhs.rep) {
        TpsRep<T>::release(rep);
        rep = rhs.rep->grab();
    }

    return *this;
}


template <class T>
Tps<T> &Tps<T>::operator=(const T &rhs) {
    TpsRep<T>::release(rep);
    rep = TpsRep<T>::zero();
    rep->data()[0] = rhs;
    return *this;
}


template <class T>
Tps<T> Tps<T>::filter(int minOrder, int maxOrder) const {
    // Compute order limits.
    maxOrder = std::min(maxOrder, getMaxOrder());
    int trcOrder = getTruncOrder();
    int variables = getVariables();

    // Construct filtered TpsRep.
    TpsRep<T> *p = TpsRep<T>::create(maxOrder, trcOrder, variables);
    int bot = getSize(minOrder - 1);
    int top = getSize(maxOrder);
    std::copy(rep->data() + bot, rep->data() + top, p->data() + bot);
    return Tps<T>(p);
}


template <class T>
Tps<T> Tps<T>::truncate(int trunc) {
    return filter(0, trunc);
}


template <class T>
const T Tps<T>::operator[](int index) const {
    return rep->data()[index];
}


template <class T>
T &Tps<T>::operator[](int index) {
    unique();
    return rep->data()[index];
}


template <class T>
const T Tps<T>::operator[](const TpsMonomial &monomial) const {
    return rep->data()[monomial.getIndex()];
}


template <class T>
T &Tps<T>::operator[](const TpsMonomial &monomial) {
    unique();
    return rep->data()[monomial.getIndex()];
}


template <class T>
const T Tps<T>::getCoefficient(int index) const {
    if(index > rep->len) {
        throw CLRangeError("Tps::getCoefficients()", "Monomial index out of range.");
    }

    return rep->data()[index];
}


template <class T>
void Tps<T>::setCoefficient(int index, const T &value) {
    int order = getOrder(index);
    if(order > rep->trcOrd) return;

    if(order > rep->maxOrd) {
        TpsRep<T> *p = TpsRep<T>::create(order, rep->trcOrd, getVariables());
        std::copy(rep->data(), rep->data() + rep->len, p->data());
        TpsRep<T>::release(rep);
        rep = p;
    } else {
        unique();
    }

    rep->data()[index] = value;
}


template <class T>
const T Tps<T>::getCoefficient(const TpsMonomial &monomial) const {
    int v1 = monomial.getVariables();
    int v2 = getVariables();

    if(v1 != v2) {
        throw SizeError("Tps::getCoefficients()",
                        "Inconsistent number of variables.");
    }

    return getCoefficient(monomial.getIndex());
}


template <class T>
void Tps<T>::setCoefficient(const TpsMonomial &monomial, const T &value) {
    int v1 = monomial.getVariables();
    int v2 = getVariables();

    if(v1 != v2) {
        throw SizeError("Tps::getCoefficients()",
                        "Inconsistent number of variables.");
    }

    setCoefficient(monomial.getIndex(), value);
}


template <class T>
Tps<T> Tps<T>::makeVariable(int nVar, int var) {
    TpsRep<T> *p = TpsRep<T>::create(1, EXACT, nVar);
    p->data()[var+1] = T(1);
    return Tps<T>(p);
}


template <class T>
Tps<T> Tps<T>::makeVarPower(int nVar, int var, int order) {
    TpsMonomial monomial(nVar);
    monomial[var] = order;
    Tps<T> z = Tps<T>(TpsRep<T>::create(order, EXACT, nVar));
    z[monomial] = T(1);
    return z;
}


template <class T>
Tps<T> Tps<T>::makeMonomial(const TpsMonomial &monomial, const T &cc) {
    int order = monomial.getOrder();
    int nVar = monomial.getVariables();
    Tps<T> z = Tps<T>(TpsRep<T>::create(order, EXACT, nVar));
    z[monomial] = cc;
    return z;
}


template <class T>
Tps<T> Tps<T>::operator+() const {
    return *this;
}


template <class T>
Tps<T> Tps<T>::operator-() const {
    TpsRep<T> *p = TpsRep<T>::create
                   (getMaxOrder(), getTruncOrder(), getVariables());
    std::transform(rep->data(), rep->data() + rep->len, p->data(),
                   std::negate<T>());
    return Tps<T>(p);
}


template <class T>
Tps<T> &Tps<T>::operator+=(const Tps<T> &rhs) {
    if(int v1 = getVariables()) {
        if(int v2 = rhs.getVariables()) {
            if(v1 != v2) {
                throw SizeError("TpsRep::operator+=()",
                                "Number of variables inconsistent.");
            }

            int trunc = std::min(getTruncOrder(), rhs.getTruncOrder());
            int xOrder = std::min(getMaxOrder(), trunc);
            int yOrder = std::min(rhs.getMaxOrder(), trunc);
            int xLength = getSize(xOrder);
            int yLength = getSize(yOrder);
            int xyLength = std::min(xLength, yLength);
            TpsRep<T> *p = TpsRep<T>::create(std::max(xOrder, yOrder), trunc, v1);
            const T *x = rep->data();
            const T *y = rhs.rep->data();
            T *z = p->data();
            std::transform(x, x + xyLength, y, z, std::plus<T>());
            std::copy(x + xyLength, x + xLength, z + xyLength);
            std::copy(y + xyLength, y + yLength, z + xyLength);
            TpsRep<T>::release(rep);
            rep = p;
        } else {
            unique();
            rep->data()[0] += rhs[0];
        }
    } else {
        *this = rhs + rep->data()[0];
    }
    return *this;
}


template <class T>
Tps<T> &Tps<T>::operator-=(const Tps<T> &rhs) {
    if(int v1 = getVariables()) {
        if(int v2 = rhs.getVariables()) {
            if(v1 != v2) {
                throw SizeError("TpsRep::operator-=()",
                                "Number of variables inconsistent.");
            }

            int trunc = std::min(getTruncOrder(), rhs.getTruncOrder());
            int xOrder = std::min(getMaxOrder(), trunc);
            int yOrder = std::min(rhs.getMaxOrder(), trunc);
            int xLength = getSize(xOrder);
            int yLength = getSize(yOrder);
            int xyLength = std::min(xLength, yLength);

            TpsRep<T> *p = TpsRep<T>::create(std::max(xOrder, yOrder), trunc, v1);
            const T *x = rep->data();
            const T *y = rhs.rep->data();
            T *z = p->data();
            std::transform(x, x + xyLength, y, z, std::minus<T>());
            std::copy(x + xyLength, x + xLength, z + xyLength);
            std::transform(y + xyLength, y + yLength, z + xyLength,
                           std::negate<T>());
            TpsRep<T>::release(rep);
            rep = p;
        } else {
            unique();
            rep->data()[0] -= rhs[0];
        }
    } else {
        *this = - rhs + rep->data()[0];
    }
    return *this;
}


template <class T>
Tps<T> &Tps<T>::operator*=(const Tps<T> &rhs) {
    return *this = multiply(rhs, truncOrder);
}


template <class T>
Tps<T> &Tps<T>::operator/=(const Tps<T> &rhs) {
    return *this = multiply(rhs.inverse(truncOrder), truncOrder);
}


template <class T>
Tps<T> &Tps<T>::operator+=(const T &rhs) {
    unique();
    rep->data()[0] += rhs;
    return *this;
}


template <class T>
Tps<T> &Tps<T>::operator-=(const T &rhs) {
    unique();
    rep->data()[0] -= rhs;
    return *this;
}


template <class T>
Tps<T> &Tps<T>::operator*=(const T &rhs) {
    unique();
    T *x = rep->data();
    std::transform(x, x + rep->len, x, std::bind(std::multiplies<T>(), std::placeholders::_1, rhs));
    return *this;
}


template <class T>
Tps<T> &Tps<T>::operator/=(const T &rhs) {
    if(rhs == T(0)) throw DivideError("Tps::operator/()");
    T *x = rep->data();
    std::transform(x, x + rep->len, x, std::bind(std::divides<T>(), std::placeholders::_1, rhs));
    return *this;
}


template <class T>
bool Tps<T>::operator==(const Tps<T> &rhs) const {
    if(int v1 = getVariables()) {
        if(int v2 = rhs.getVariables()) {
            if(v1 == v2) {
                int trunc = std::min(getTruncOrder(), rhs.getTruncOrder());
                int xOrder = std::min(getMaxOrder(), trunc);
                int yOrder = std::min(rhs.getMaxOrder(), trunc);
                int xLength  = getSize(xOrder);
                int yLength  = getSize(yOrder);
                int xyLength = getSize(std::min(xOrder, yOrder));
                const T *x = rep->data();
                const T *y = rhs.rep->data();

                for(int i = 0; i < xyLength; i++) {
                    if(x[i] != y[i]) return false;
                }

                for(int i = xyLength; i < xLength; i++) {
                    if(x[i] != T(0)) return false;
                }

                for(int i = xyLength; i < yLength; i++) {
                    if(y[i] != T(0)) return false;
                }

                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    } else {
        if(rhs.getVariables()) {
            return false;
        } else {
            return rep->data()[0] == rhs.rep->data()[0];
        }
    }
}


template <class T>
bool Tps<T>::operator==(const T &rhs) const {
    const T *x = rep->data();

    if(x[0] != rhs) return false;

    for(int i = 1; i < getSize(); ++i) {
        if(x[i] != T(0)) return false;
    }

    return true;
}


template <class T>
bool Tps<T>::operator!=(const Tps<T> &rhs) const {
    return !(*this == rhs);
}


template <class T>
bool Tps<T>::operator!=(const T &rhs) const {
    return !(*this == rhs);
}


template <class T>
Tps<T> Tps<T>::substitute(const Matrix<T> &M) const {
    if(getVariables()) {
        int v1 = getVariables();
        int v2 = M.nrows();
        if(v1 != v2) {
            throw SizeError("Tps::substitute()", "Matrix not consistent with Tps.");
        }
        int nRow = M.nrows();
        int nCol = M.ncols();

        // Define the nRow linear transformations.
        Array1D< Tps<T> > y(nRow);

        for(int i = 0; i < nRow; ++i) {
            y[i] = Tps<T>(1, nCol);
            for(int j = 0; j < nCol; ++j) y[i][j+1] = M[i][j];
        }

        // Evaluate the substitution.
        const T *x = rep->data();
        Tps<T> z(x[0]);

        if(int maxOrd = getMaxOrder()) {
            const Array1D<TpsSubstitution> &table = rep->help->getSubTable();
            Array1D< Tps<T> > product(maxOrd + 1);
            product[0] = Tps<T>(T(1));

            for(int next = 1; next < table.size();) {
                const TpsSubstitution &s = table[next];
                product[s.order] = product[s.order-1] * y[s.variable];
                z += x[s.index] * product[s.order];
                next = (s.order < maxOrd) ? next + 1 : s.skip;
            }
        }

        return z;
    } else {
        return *this;
    }
}


template <class T>
Tps<T> Tps<T>::substitute(const VpsMap<T> &rhs) const {
    int v1 = getVariables();
    int v2 = rhs.getDimension();
    if(v1 != v2) {
        throw SizeError("Tps::substitute()", "VpsMap is inconsistent with Tps.");
    }

    const T *x = rep->data();
    Tps<T> z(x[0]);

    if(int maxOrd = getMaxOrder()) {
        const Array1D<TpsSubstitution> &table = rep->help->getSubTable();
        Array1D< Tps<T> > product(maxOrd + 1);
        product[0] = Tps<T>(T(1));
        int trunc = getTruncOrder();

        for(int next = 1; next < table.size();) {
            const TpsSubstitution &s = table[next];
            product[s.order] = product[s.order-1].multiply(rhs[s.variable], trunc);
            z += x[s.index] * product[s.order];
            next = (s.order < maxOrd) ? next + 1 : s.skip;
        }
    }

    return z;
}


template <class T>
T Tps<T>::evaluate(const Vector<T> &rhs) const {
    int v1 = getVariables();
    int v2 = rhs.size();
    if(v1 != v2) {
        throw SizeError("Tps::evaluate()", "Vector is inconsistent with Tps.");
    }

    const T *x = rep->data();
    T z = x[0];

    if(int maxOrd = getMaxOrder()) {
        const Array1D<TpsSubstitution> &table = rep->help->getSubTable();
        Array1D<T> product(maxOrd + 1);
        product[0] = T(1);

        for(int next = 1; next < table.size();) {
            const TpsSubstitution &s = table[next];
            product[s.order] = product[s.order-1] * rhs[s.variable];
            z += x[s.index] * product[s.order];
            next = (s.order < maxOrd) ? next + 1 : s.skip;
        }
    }

    return z;
}


template <class T>
void Tps<T>::clear() {
    TpsRep<T>::release(rep);
    rep = TpsRep<T>::zero();
}


template <class T>
std::istream &Tps<T>::get(std::istream &is) {
    is.flags(std::ios::skipws);
    char head[4];
    is.get(head, 4);
    if(strcmp(head, "Tps") != 0) {
        throw FormatError("Tps::get()", "Flag word \"Tps\" missing.");
    }

    int maxOrder, truncOrder, nVar;
    is >> maxOrder >> truncOrder >> nVar;
    Tps<T> z(TpsRep<T>::create(maxOrder, truncOrder, nVar));
    T coeff;

    if(nVar <= 0  ||  truncOrder == 0) {
        is >> coeff;
        z[0] = coeff;

        if(coeff != T(0)) {
            z[0] = coeff;
            is >> coeff;
        }
    } else {
        TpsMonomial monomial(nVar);
        maxOrder = 0;
        bool done = false;
        bool fail = false;

        while(true) {
            is >> coeff;
            fail = is.fail();

            int order = 0;
            for(int var = 0; var < nVar; var++) {
                int p;
                is >> p;
                fail |= is.fail();
                if(p < 0) done = true;
                monomial[var] = p;
                order += monomial[var];
            }

            if(done) break;
            if(fail) throw FormatError("Tps::get()", "File read error");
            int index = monomial.getIndex();

            if(coeff != T(0)) {
                maxOrder = order;
            } else if(index == 0) {
                break;
            }

            z[index] = coeff;
        }

        z.setMaxOrder(maxOrder);
        *this = z;
    }

    return is;
}


template <class T>
std::ostream &Tps<T>::put(std::ostream &os) const {
    std::streamsize old_prec = os.precision(14);
    os.setf(std::ios::scientific, std::ios::floatfield);

    int nVar = getVariables();
    os << "Tps " << getMaxOrder() << ' ' << getTruncOrder() << ' '
       << nVar << std::endl;

    if(nVar == 0) {
        os << std::setw(24) << rep->data()[0] << std::endl;
    } else {
        for(int i = 0; i < getSize(); ++i) {
            if(rep->data()[i] != T(0)) {
                os << std::setw(24) << rep->data()[i];

                for(int var = 0; var < nVar; var++) {
                    os << std::setw(3) << getExponents(i)[var];
                }

                os << std::endl;
            }
        }

        os << std::setw(24) << T(0);

        for(int var = 0; var < nVar; var++) {
            os << std::setw(3) << (-1);
        }
    }

    os << std::endl;

    os.precision(old_prec);
    os.setf(std::ios::fixed, std::ios::floatfield);
    return os;
}


template <class T>
Tps<T> Tps<T>::multiply(const Tps<T> &rhs, int trunc) const {
    int v1 = getVariables();
    int v2 = rhs.getVariables();
    if(v1) {
        if(v2) {
            if(v1 != v2) {
                throw SizeError("TpsRep::multiply()",
                                "Number of variables inconsistent.");
            }

            if(getTruncOrder() != EXACT) {
                int cut = getTruncOrder();
                if(rhs[0] == 0.0) ++cut;
                trunc = std::min(trunc, cut);
            }

            if(rhs.getTruncOrder() != EXACT) {
                int cut = rhs.getTruncOrder();
                if((*this)[0] == 0.0) ++cut;
                trunc = std::min(trunc, cut);
            }

            int maxOrder = std::min(getMaxOrder() + rhs.getMaxOrder(), trunc);

            TpsRep<T> *p = TpsRep<T>::create(maxOrder, trunc, v1);
            const T *x = rep->data();
            T *z = p->data();
            int yBot = 0;
            int yHig = std::min(rhs.getMaxOrder(), trunc);

            for(int yOrd = 0; yOrd <= yHig; yOrd++) {
                int xOrd = std::min(getMaxOrder(), trunc - yOrd);
                int xTop = getSize(xOrd);
                int yTop = getSize(yOrd);

                for(int yInd = yBot; yInd < yTop; yInd++) {
                    T y = rhs.rep->data()[yInd];
                    if(y != T(0)) {
                        const int *prod = rep->help->getProductArray(yInd);

                        for(int xInd = 0; xInd < xTop; xInd++) {
                            z[prod[xInd]] += x[xInd] * y;
                        }
                    }
                }

                yBot = yTop;
            }

            return Tps<T>(p);
        } else {
            Tps<T> result(TpsRep<T>::create
                          (getMaxOrder(), getTruncOrder(), v1));
            const T *x = rep->data();
            const T y  = rhs.rep->data()[0];
            T *z = result.rep->data();
            std::transform(x, x + rep->len, z,
                           std::bind(std::multiplies<T>(), std::placeholders::_1, y));
            return result;
        }
    } else {
        Tps<T> result(TpsRep<T>::create
                      (rhs.getMaxOrder(), rhs.getTruncOrder(), v2));
        const T x  = rep->data()[0];
        const T *y = rhs.rep->data();
        T *z = result.rep->data();
        std::transform(y, y + rep->len, z,
                       std::bind(std::multiplies<T>(), std::placeholders::_1, x));
        return result;
    }
}


template <class T>
Tps<T> Tps<T>::inverse(int trunc) const {
    T aZero = rep->data()[0];
    if(aZero == T(0)) throw DivideError("Tps::inverse()");

    if(isConstant()) {
        return Tps<T>(T(1) / aZero);
    } else {
        int cut = std::min(trunc, getTruncOrder());
        T *series = new T[cut+1];
        series[0] = T(1) / aZero;

        for(int i = 1; i <= cut; i++) {
            series[i] = - series[i-1] / aZero;
        }

        Tps<T> z = Taylor(series, cut);
        delete [] series;
        return z;
    }
}


template <class T>
Tps<T> Tps<T>::derivative(int var) const {
    if(getVariables()  &&  getMaxOrder() > 0) {
        int maxOrder = getMaxOrder() - 1;
        int trcOrder = getTruncOrder();

        TpsRep<T> *p = TpsRep<T>::create(maxOrder, trcOrder, getVariables());
        const T *x = rep->data();
        T *z = p->data();

        const int *product = rep->help->getProductArray(var + 1);

        for(int i = getSize(maxOrder); i-- > 0;) {
            int k = product[i];
            z[i] = x[k] * double(getExponents(k)[var]);
        }

        Tps<T> result(p);
        return result;
    } else {
        Tps<T> result(TpsRep<T>::zero());
        return result;
    }
}


template <class T>
Tps<T> Tps<T>::integral(int var) const {
    if(getVariables() == 0) {
        throw LogicalError("TpsRep::integral()", "Cannot integrate a constant.");
    }

    int trcO = std::min(rep->trcOrd + 1, truncOrder);
    int maxO = std::min(rep->maxOrd + 1, trcO);
    TpsRep<T> *p = TpsRep<T>::create(maxO, trcO, getVariables());

    const T *x = rep->data();
    T *z = p->data();
    const int *product = rep->help->getProductArray(var + 1);

    for(int i = getSize(maxO - 1); i-- > 0;) {
        int k = product[i];
        z[k] = x[i] / double(getExponents(k)[var]);
    }

    return Tps<T>(p);
}


template <class T>
Tps<T> Tps<T>::multiplyVariable(int var) const {
    if(getVariables() == 0) {
        throw LogicalError("TpsRep::multiplyVariable()",
                           "Cannot multiply a constant by a numbered variable.");
    }

    int trcO = std::min(rep->trcOrd + 1, truncOrder);
    int maxO = std::min(rep->maxOrd + 1, trcO);
    TpsRep<T> *p = TpsRep<T>::create(maxO, trcO, getVariables());

    const T *x = rep->data();
    T *z = p->data();
    const int *product = rep->help->getProductArray(var + 1);

    for(int i = getSize(maxO - 1); i-- > 0;) {
        z[product[i]] = x[i];
    }

    return Tps<T>(p);
}


template <class T>
Tps<T> Tps<T>::scaleMonomials(const Tps<T> &rhs) const {
    int v1 = getVariables();
    int v2 = rhs.getVariables();

    if(v1 != v2) {
        throw SizeError("TpsRep::scaleMonomials()",
                        "Number of variables inconsistent.");
    }

    int order = std::min(getMaxOrder(), rhs.getMaxOrder());
    int trunc = std::min(getTruncOrder(), rhs.getTruncOrder());

    TpsRep<T> *p = TpsRep<T>::create(std::min(order, trunc), trunc, v1);
    const T *x = rep->data();
    const T *y = rhs.rep->data();
    T *z = p->data();
    std::transform(x, x + p->len, y, z, std::multiplies<T>());
    return Tps<T>(p);
}


template <class T>
Tps<T> Tps<T>::Taylor(const T series[], int order) const {
    if(isConstant()) {
        return Tps<T>(series[0]);
    } else {
        Tps<T> x(*this);
        x[0] = T(0);
        Tps<T> z(series[order]);
        for(int maxOrder = 1; maxOrder <= order; maxOrder++) {
            z = x.multiply(z, maxOrder);
            z[0] = series[order-maxOrder];
        }
        return z;
    }
}


template <class T>
int Tps<T>::getMaxOrder() const {
    return rep->maxOrd;
}


template <class T>
int Tps<T>::getTruncOrder() const {
    return std::min(rep->trcOrd, truncOrder);
}


template <class T>
int Tps<T>::getVariables() const {
    return (rep->help != 0) ? rep->help->getVariables() : 0;
}


template <class T>
int Tps<T>::getSize() const {
    return rep->len;
}


template <class T>
bool Tps<T>::isConstant() const {
    return rep->help == 0;
}


template <class T>
int Tps<T>::getGlobalTruncOrder() {
    return truncOrder;
}


template <class T>
void Tps<T>::setGlobalTruncOrder(int order) {
    truncOrder = order;
}


template <class T>
const TpsMonomial &Tps<T>::getExponents(int index) const {
    if(rep->help == 0) {
        throw LogicalError("Tps::getExponents()",
                           "Cannot get exponents of a constant.");
    }

    return rep->help->getExponents(index);
}


template <class T>
int Tps<T>::getOrder(int index) const {
    return rep->help ? rep->help->getOrder(index) : 0;
}


template <class T>
int Tps<T>::getSize(int order) const {
    return rep->help ? rep->help->getSize(order) : 1;
}


template <class T> inline
Tps<T>::Tps(TpsRep<T> *d): rep(d)
{}

#endif // CLASSIC_Tps_CC
