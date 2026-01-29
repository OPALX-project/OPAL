#ifndef CLASSIC_TpsMonomial_HH
#define CLASSIC_TpsMonomial_HH

// ------------------------------------------------------------------------
// $RCSfile: TpsMonomial.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TpsMonomial
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2002/03/25 20:44:15 $
// $Author: dbruhwil $
//
// ------------------------------------------------------------------------

#include "Algebra/Array1D.h"


// Class TpsMonomial
// ------------------------------------------------------------------------
/// Exponent array for Tps<T>.
//  A representation of the exponents for a Tps monomial.

class TpsMonomial {

public:

    /// Constructor.
    //  Array of zeros with length [b]nVar[/b], defines [b]nVar[/b] variables
    //  with zero exponent each.
    TpsMonomial(int nVar);

    /// Constructor.
    //  Array with length [b]nVar[/b], define variable with index [b]var[/b].
    TpsMonomial(int nVar, int var);

    TpsMonomial();
    TpsMonomial(const TpsMonomial &);
    ~TpsMonomial();
    const TpsMonomial &operator=(const TpsMonomial &);

    /// Get exponent.
    //  Return a reference to the exponent of variable [b]index[/b].
    int &operator[](int index);

    /// Get exponent.
    //  Return value of the exponent of variable [b]index[/b].
    int operator[](int index) const;

    /// Product.
    //  Return the exponent set for the product of the monomials denoted by
    //  [b]this[/b] and [b]rhs[/b], i.e. the sum of the exponents for each
    //  variable.
    TpsMonomial operator*(const TpsMonomial &rhs) const;

    /// Convert.
    //  Convert the monomial to the corresponding Giorgilli index.
    int getIndex() const;

    /// Get order.
    //  Compute the monomial's order, i.e. the sum of all exponents.
    int getOrder() const;

    /// Get variables.
    //  Return the monomial's number of variables, i.e. the number of exponents.
    int getVariables() const;

private:

    // The monomial's exponents.
    Array1D<int> array;
};

#endif // CLASSIC_TpsMonomial_HH
