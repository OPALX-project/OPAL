// ------------------------------------------------------------------------
// $RCSfile: Vps.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Vps
//   Vector truncated power series in n variables.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/Vps.h"
#include "Algebra/Matrix.h"
#include "Algebra/LUMatrix.h"
#include "Algebra/Vector.h"
#include "Algebra/Tps.h"
#include "Algebra/TpsMonomial.h"
#include "Utilities/DivideError.h"
#include "Utilities/FormatError.h"
#include "Utilities/SizeError.h"
#include <algorithm>
#include <iostream>

#include <cstring>


// Template class Vps<T>
// ------------------------------------------------------------------------

template <class T>
Vps<T>::Vps():
    data(), variables(0)
{}


template <class T>
Vps<T>::Vps(int nDim, int nVar):
    data(nDim, Tps<T>(0, nVar)), variables(nVar)
{}


template <class T>
Vps<T>::Vps(const Vps<T> &rhs):
    data(rhs.data), variables(rhs.variables) {
    // No call to check, assume rhs is correct.
}


template <class T>
Vps<T>::Vps(const Matrix<T> &x):
    data(x.nrows()), variables(x.ncols()) {
    for(int i = 0; i < x.nrows(); i++) {
        data[i] = Tps<T>(1, x.ncols());

        for(int j = 0; j < x.ncols(); j++) {
            data[i][j+1] = x(i, j);
        }
    }
}


template <class T>
Vps<T>::Vps(const Vector<T> &x):
    data(x.size()), variables(0) {
    for(int i = 0; i < getDimension(); i++) data[i] = Tps<T>(x[i]);
}


template <class T>
Vps<T>::~Vps()
{}


template <class T>
Vps<T> &Vps<T>::operator=(const Vps<T> &y) {
    data = y.data;
    variables = y.variables;
    return *this;
}


template <class T>
const Tps<T> &Vps<T>::getComponent(int index) const {
    return data(index);
}


template <class T>
void Vps<T>::setComponent(int index, const Tps<T> &value) {
    check();
    if(value.getVariables() != variables) {
        throw SizeError("Vps::setComponent()", "Index out of range.");
    }
    data(index) = value;
}


template <class T>
void Vps<T>::check() const {
    variables = 0;

    for(const Tps<T> *p = data.begin(); p != data.end(); ++p) {
        if(int var = p->getVariables()) {
            if(variables != 0  &&  var != variables) {
                throw SizeError("Vps::check()",
                                "Vps has inconsistent number of variables.");
            }
            variables = var;
        }
    }
}


template <class T>
Vps<T> Vps<T>::operator+() const {
    return *this;
}


template <class T>
Vps<T> Vps<T>::operator-() const {
    Vps<T> z(getDimension(), variables);
    for(int i = 0; i < getDimension(); i++) z[i] = - data[i];
    return z;
}


template <class T>
Vps<T> &Vps<T>::operator*=(const Tps<T> &y) {
    for(int i = 0; i < getDimension(); i++) data[i] *= y;
    return *this;
}


template <class T>
Vps<T> &Vps<T>::operator/=(const Tps<T> &y) {
    Tps<T> t = y.inverse(Tps<T>::getGlobalTruncOrder());
    for(int i = 0; i < getDimension(); i++) data[i] *= t;
    return *this;
}


template <class T>
Vps<T> &Vps<T>::operator*=(const T &y) {
    for(int i = 0; i < getDimension(); i++) data[i] *= y;
    return *this;
}


template <class T>
Vps<T> &Vps<T>::operator/=(const T &y) {
    if(y == T(0)) throw DivideError("Vps::operator/()");
    T t = T(1) / y;
    for(int i = 0; i < getDimension(); i++) data[i] *= t;
    return *this;
}


template <class T>
Vps<T>& Vps<T>::operator+=(const Vps<T> &y) {
    if(getDimension() != y.getDimension()) {
        throw SizeError("Vps::operator+=()", "Inconsistent dimensions.");
    }

    for(int i = 0; i < getDimension(); i++) data[i] += y.data[i];
    return *this;
}


template <class T>
Vps<T>& Vps<T>::operator-=(const Vps<T> &y) {
    if(getDimension() != y.getDimension()) {
        throw SizeError("Vps::operator-=()", "Inconsistent dimensions.");
    }

    for(int i = 0; i < getDimension(); i++) data[i] -= y.data[i];
    return *this;
}


template <class T>
Vps<T>& Vps<T>::operator+=(const Vector<T> &y) {
    if(getDimension() != y.size()) {
        throw SizeError("Vps::operator+=()", "Inconsistent dimensions.");
    }

    for(int i = 0; i < getDimension(); i++) data[i] += y[i];
    return *this;
}


template <class T>
Vps<T>& Vps<T>::operator-=(const Vector<T> &y) {
    if(getDimension() != y.size()) {
        throw SizeError("Vps::operator-=()", "Inconsistent dimensions.");
    }

    for(int i = 0; i < getDimension(); i++) data[i] -= y[i];
    return *this;
}


template <class T>
std::istream &Vps<T>::get(std::istream &is) {
    char head[4];
    is.flags(std::ios::skipws);
    is.get(head, 4);

    if(std::strcmp(head, "Vps") != 0) {
        throw FormatError("Vps::get()", "Flag word \"Vps\" missing.");
    }

    int nDim;
    is >> nDim;

    if(nDim == 0) {
        throw FormatError("Vps::get()", "Zero Vps dimension");
    }

    Vps<T> temp(nDim);
    for(int i = 0; i < nDim; i++) is >> data[i];
    temp.check();
    *this = temp;
    return is;
}


template <class T>
std::ostream &Vps<T>::put(std::ostream &os) const {
    int nDim = getDimension();
    os << "Vps " << nDim << std::endl;
    for(int i = 0; i < nDim; i++) os << data[i];
    return os;
}


template <class T>
int Vps<T>::getDimension() const {
    return data.size();
}


template <class T>
int Vps<T>::getTopOrder() const {
    const Tps<T> *p = data.begin();
    int topOrder = p->getMaxOrder();

    while(++p < data.end()) {
        topOrder = std::max(topOrder, p->getMaxOrder());
    }

    return topOrder;
}


template <class T>
int Vps<T>::getTruncOrder() const {
    const Tps<T> *p = data.begin();
    int truncOrder = p->getTruncOrder();

    while(++p < data.end()) {
        truncOrder = std::min(truncOrder, p->getTruncOrder());
    }

    return truncOrder;
}


template <class T>
int Vps<T>::getVariables() const {
    check();
    return variables;
}


template <class T>
Vps<T> Vps<T>::filter(int lowOrder, int highOrder) const {
    check();
    Vps<T> z(getDimension(), variables);

    for(int i = 0; i < getDimension(); i++) {
        z.data[i] = data[i].filter(lowOrder, highOrder);
    }

    return z;
}


template <class T>
Vps<T> Vps<T>::truncate(int trunc) {
    check();
    Vps<T> z(getDimension(), variables);

    for(int i = 0; i < getDimension(); i++) {
        z.data[i] = data[i].truncate(trunc);
    }

    return *this;
}
