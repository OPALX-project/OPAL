// ------------------------------------------------------------------------
// $RCSfile: SFunction.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: SFunction
//   Object handling the position funcitions SI(), SC(), and SO().
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:42 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Expressions/SFunction.h"
#include "Utilities/OpalException.h"


// The pointer to the current s-dependent function.
// ------------------------------------------------------------------------

const SFunction *SFunction::sfun = 0;


// Class SFunction
// ------------------------------------------------------------------------

SFunction::SFunction() {
    reset();
    sfun = this;
}


SFunction::~SFunction() {
    sfun = 0;
}


double SFunction::arcIn() {
    if(sfun) {
        return sfun->position(1.0);
    } else {
        throw OpalException("arcIn()",
                            "The use of function \"SI()\" "
                            "is not valid in this context.");
    }
}


double SFunction::arcCtr() {
    if(sfun) {
        return sfun->position(0.5);
    } else {
        throw OpalException("arcCtr()",
                            "The use of function \"SC()\" "
                            "is not valid in this context.");
    }
}


double SFunction::arcOut() {
    if(sfun) {
        return sfun->position(0.0);
    } else {
        throw OpalException("arcOut()",
                            "The use of function \"SO()\" "
                            "is not valid in this context.");
    }
}


void SFunction::reset() {
    elementLength = 0.0;
    exitArc = 0.0;
}


void SFunction::update(double length) {
    elementLength = length;
    exitArc += elementLength;
}


double SFunction::position(double flag) const {
    // Backtrack to desired position.
    return exitArc - flag * elementLength;
}
