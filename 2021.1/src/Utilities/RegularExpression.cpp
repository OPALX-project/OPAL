// ------------------------------------------------------------------------
// $RCSfile: RegularExpression.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RegularExpression.
//   This class provides a simple interface to the POSIX-compliant
//   regcomp/regexec package.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Utilities/RegularExpression.h"
#include "Utilities/OpalException.h"
#include <cstdlib>
#include <regex.h>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <new>
#include <cstring>
#include <functional>

// Class RegularExpression::Expression
// ------------------------------------------------------------------------

class RegularExpression::Expression: public regex_t {};


// Class RegularExpression
// ------------------------------------------------------------------------

RegularExpression::RegularExpression(const std::string &pattern, bool ignore):
    patt(pattern), caseIgnore(ignore) {
    init();
}


RegularExpression::RegularExpression(const RegularExpression &rhs):
    patt(rhs.patt), caseIgnore(rhs.caseIgnore) {
    init();
}


RegularExpression::~RegularExpression() {
    ::regfree(expr);
}


bool RegularExpression::match(const std::string &s) const {
    int rc = state;

    if(rc == 0) {
        int nmatch = 0;
        regmatch_t *pmatch = 0;
        int flags = 0;
        rc = regexec(expr, s.data(), nmatch, pmatch, flags);
        if(rc == 0) {
            return true;
        } else if(rc == REG_NOMATCH) {
            return false;
        }
    }

    char buffer[80];
    regerror(rc, expr, buffer, 80);
    throw OpalException("RegularExpression::match()", buffer);
}


bool RegularExpression::OK() const {
    return state == 0  &&  expr != 0;
}


void RegularExpression::init() {
    expr = (Expression *) malloc(sizeof(regex_t));
    int flags = REG_NOSUB;
    if(caseIgnore) flags |= REG_ICASE;
    state = regcomp(expr, patt.c_str(), flags);
}
