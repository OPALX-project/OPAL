// ------------------------------------------------------------------------
// $RCSfile: TpsData.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TpsData
//   Bookkeeper class for Tps class.
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/TpsData.h"
#include "Algebra/TpsMonomial.h"


// Class TpsData
// ------------------------------------------------------------------------

namespace {
    const int TpsData_SIZE = 100;
    TpsData theirTpsData[TpsData_SIZE];
};


TpsData::TpsData() {
    clear();
}


TpsData::~TpsData()
{}


TpsData *TpsData::getTpsData(int nOrd, int nVar) {
    if(nVar == 0) {
        return 0;
    } else {
        theirTpsData[nVar].build(nOrd, nVar);
        return &theirTpsData[nVar];
    }
}


int TpsData::indexMonomial(const TpsMonomial &pows) const {
    int index = 0;
    int order = 0;

    for(int var = variables; var-- > 0;) {
        order += pows[var];
        index += binom[var][order];
    }

    return index;
}


void TpsData::build(int nOrd, int nVar) {
    // Number of variables.
    variables = nVar;

    // Maximum order so far.
    if(topOrder < nOrd) {
        topOrder = nOrd;

        // Define binomial coefficients.
        int var, ord;
        binom = BinomialTable(nVar, topOrder + 2);

        for(var = 0; var < nVar; var++) {
            binom[var][0] = 0;
        }

        for(ord = 1; ord <= topOrder + 1; ord++) {
            binom[nVar-1][ord] = ord;

            for(var = nVar - 1; var-- > 0;) {
                binom[var][ord] = binom[var][ord-1] + binom[var+1][ord];
            }
        }

        // Set up table of exponents up to topOrder
        // (algorithm due to Liam Healey (Marylie 3.0)).
        maxSize = getSize(topOrder);
        {
            TpsMonomial power(nVar);
            expon = ExponentTable(maxSize, power);

            for(int index = 1; index < maxSize; index++) {
                int carry = power[nVar-1];
                power[nVar-1] = 0;
                int lastnz = int(-1);

                for(int j = 0; j < nVar - 1; j++) {
                    if(power[j] > 0) lastnz = j;
                }

                if(lastnz == int(-1)) {
                    power[lastnz+1] = power[lastnz+1] + 1 + carry;
                } else {
                    power[lastnz]--;
                    power[lastnz+1] = power[lastnz+1] + 1 + carry;
                }

                expon[index] = power;
            }
        }

        // Set up table of products for total order up to topOrder.
        {
            prod = ProductTable(maxSize);
            TpsMonomial power(nVar);

            for(int xord = 0; xord <= topOrder; xord++) {
                int yord = topOrder - xord;
                int ysize = getSize(yord);
                for(int i = getSize(xord - 1); i < getSize(xord); i++) {
                    prod[i] = ProductRow(ysize);

                    for(int j = 0; j < ysize; j++) {
                        power = expon[i] * expon[j];
                        prod[i][j] = indexMonomial(power);
                    }
                }
            }
        }

        {
            // Fill the substitution table.
            subTable = Array1D<TpsSubstitution>(maxSize);
            TpsMonomial power(nVar);
            int next = 1;
            fillSubst(0, 1, power, next);
        }
    }
}


void TpsData::clear() {
    variables = 0;
    topOrder = 0;
}


void TpsData::fillSubst
(int var, int order, TpsMonomial &pow, int &next) {
    for(int v = var; v < variables; v++) {
        TpsSubstitution &s = subTable[next++];
        pow[v]++;
        s.index = indexMonomial(pow);
        s.order = order;
        s.variable = v;
        if(order < topOrder) fillSubst(v, order + 1, pow, next);
        s.skip = next;
        pow[v]--;
    }
}
