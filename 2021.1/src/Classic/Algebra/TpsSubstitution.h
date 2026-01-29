#ifndef CLASSIC_TpsSubstitution_H
#define CLASSIC_TpsSubstitution_H

// ------------------------------------------------------------------------
// $RCSfile: TpsSubstitution.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Struct: TpsSubstitution
//
// ------------------------------------------------------------------------
// Class category: Algebra
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------


// Struct TpsSubstitution
// -------------------------------------------------------------------------
/// Substitution for Tps<T>.
//  This is an internal bookkeeping class, used for substitution into
//  a truncated power series.  It should not be used by user programs.


struct TpsSubstitution {

    TpsSubstitution();
    TpsSubstitution(const TpsSubstitution &);
    ~TpsSubstitution();
    void operator=(const TpsSubstitution &);

    int index;
    int order;
    int variable;
    int skip;
};


inline TpsSubstitution::TpsSubstitution():
    index(0), order(0), variable(0), skip(0)
{}


inline TpsSubstitution::TpsSubstitution(const TpsSubstitution &rhs):
    index(rhs.index), order(rhs.order), variable(rhs.variable), skip(rhs.skip)
{}


inline TpsSubstitution::~TpsSubstitution()
{}


inline void TpsSubstitution::operator=(const TpsSubstitution &rhs) {
    index = rhs.index;
    order = rhs.order;
    variable = rhs.variable;
    skip = rhs.skip;
}

#endif // CLASSIC_TpsSubstitution_H
