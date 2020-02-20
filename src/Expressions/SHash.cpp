// ------------------------------------------------------------------------
// $RCSfile: SHash.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SHash
//   Representation of a "#" value in a ATable array expression.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Expressions/SHash.h"
#include "Expressions/ATable.h"
#include "Expressions/SConstant.h"
#include <iostream>


// Class SHash
// ------------------------------------------------------------------------

namespace Expressions {

    SHash::SHash(const SHash &rhs):
        Scalar<double>(rhs),
        itsTable(rhs.itsTable)
    {}


    SHash::SHash(const ATable &table):
        itsTable(table)
    {}


    SHash::~SHash()
    {}


    Scalar<double> *SHash::clone() const {
        return new SConstant<double>(itsTable.getHash());
    }


    double SHash::evaluate() const {
        return itsTable.getHash();
    }


    void SHash::print(std::ostream &os, int /*precedence*/) const {
        os << '#';
    }

}
