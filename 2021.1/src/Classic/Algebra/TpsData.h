#ifndef CLASSIC_TpsData_H
#define CLASSIC_TpsData_H

// ------------------------------------------------------------------------
// $RCSfile: TpsData.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TpsData
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"
#include "Algebra/Array2D.h"
#include "Algebra/TpsMonomial.h"
#include "Algebra/TpsSubstitution.h"


// Class TpsData
// ------------------------------------------------------------------------
/// Bookkeeping class for Tps<T>.
//  Internal utility class for Tps class.
//  Not to be used directly in a program.

class TpsData {

public:

    TpsData();
    ~TpsData();

    // Return the TpsData structure for "nVar" variables and order "nOrd".
    static TpsData *getTpsData(int nOrd, int nVar);

    // Power array for monomial with index "index"
    const TpsMonomial &getExponents(int index) const;

    // Order of monomial with index "index".
    int getOrder(int index) const;

    // Index array for products of monomial "index".
    const int *getProductArray(int index) const;

    // Number of terms in a Tps object of order "order"
    int getSize(int order) const;

    // Number of terms in a Tps object of order "order"
    const Array1D<TpsSubstitution> &getSubTable() const;

    // Number of variables
    int getVariables() const;

    // Index for a monomial from its exponents
    int indexMonomial(const TpsMonomial &) const;

private:

    // Not implemented.
    TpsData(const TpsData &);
    void operator=(const TpsData &);

    // TpsData builder.
    void build(int max, int nVar);

    // Initialise TpsData to empty.
    void clear();

    // Initialise Substitution tables.
    void fillSubst(int var, int order, TpsMonomial &pow, int &next);

    // Number of variables.
    int variables;

    // Maximum order to which the data are expanded.
    int topOrder;

    // Maximum Tps size.
    int maxSize;

    // Exponents: expon[index] = set of exponents in monomial "index".
    typedef Array1D<TpsMonomial> ExponentTable;
    ExponentTable expon;

    // Indexing information, used for Tps multiplication.
    typedef Array1D<int> ProductRow;
    typedef Array1D<ProductRow> ProductTable;
    ProductTable prod;

    // Binomial coefficients: binom[i-1][k]=(i+k) ! / i! / k!
    // used in the indexMonomial() and getSize(int) methods.
    typedef Array2D<int> BinomialTable;
    BinomialTable binom;

    // The substitution table.
    Array1D<TpsSubstitution> subTable;
};


// Class TpsData, inline methods.
// ------------------------------------------------------------------------

inline const TpsMonomial &TpsData::getExponents(int index) const
{ return expon[index]; }


inline const int *TpsData::getProductArray(int index) const
{ return prod[index].begin(); }


inline int TpsData::getOrder(int index) const
{ return expon[index].getOrder(); }


inline int TpsData::getSize(int order) const
{ return topOrder == 0 ? 1 : binom[0][order+1]; }


inline const Array1D<TpsSubstitution> &TpsData::getSubTable() const
{ return subTable; }


inline int TpsData::getVariables() const
{ return variables; }

#endif // CLASSIC_TpsData_H
