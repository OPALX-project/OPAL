// ------------------------------------------------------------------------
// $RCSfile: TpsMonomial.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TpsMonomial
//   Index set for a Tps monomial.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/TpsMonomial.h"
#include "Algebra/TpsData.h"


// Class TpsMonomial
// ------------------------------------------------------------------------

TpsMonomial::TpsMonomial():
    array()
{}


TpsMonomial::TpsMonomial(const TpsMonomial &rhs):
    array(rhs.array)
{}


TpsMonomial::TpsMonomial(int nVar):
    array(nVar, int(0))
{}


TpsMonomial::TpsMonomial(int nVar, int var):
    array(nVar, int(0)) {
    array[var] = 1;
}


TpsMonomial::~TpsMonomial()
{}


const TpsMonomial &TpsMonomial::operator=(const TpsMonomial &rhs) {
    array = rhs.array;
    return *this;
}


int &TpsMonomial::operator[](int index) {
    return array[index];
}


int TpsMonomial::operator[](int index) const {
    return array[index];
}


TpsMonomial TpsMonomial::operator*(const TpsMonomial &rhs) const {
    int n = rhs.getVariables();
    TpsMonomial z(*this);
    for(int i = 0; i < n; i++) z.array[i] += rhs.array[i];
    return z;
}


int TpsMonomial::getIndex() const {
    TpsData *data = TpsData::getTpsData(getOrder(), array.size());
    return data->indexMonomial(*this);
}


int TpsMonomial::getOrder() const {
    int order = 0;

    for(int i = 0; i < array.size(); i++) {
        order += array[i];
    }

    return order;
}


int TpsMonomial::getVariables() const {
    return array.size();
}
