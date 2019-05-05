// ------------------------------------------------------------------------
// $RCSfile: Match.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class Match
//   This class encapsulates all data for matching which do not depend
//   on the matching method.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/29 10:42:09 $
// $Author: opal $
//
// ------------------------------------------------------------------------

#include "Match/Match.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/OpalData.h"
#include "Algebra/Vector.h"
#include "Match/AbstractFun.h"
#include "Match/AbstractVar.h"
#include <iomanip>
#include <iostream>

using std::cout;


// Class Match
// ------------------------------------------------------------------------

Match *Match::block = 0;


Match::Match(): parser(), constraintCount(0), callCount(0), printLevel(1)
{}


Match::~Match() {
    while(! theVariables.empty()) {
        delete theVariables.back();
        theVariables.pop_back();
    }

    while(! theFunctions.empty()) {
        delete theFunctions.back();
        theFunctions.pop_back();
    }
}


void Match::addVariable(AbstractVar *var) {
    theVariables.push_back(var);
}


void Match::deleteVariable(const std::string &name) {
    for(VarList::iterator var = theVariables.begin();
        var != theVariables.end();) { // for loop without increment
        if((*var)->getName() == name) {
            var = theVariables.erase(var); // new value: one after the erased element
        } else {
            ++var;
        }
    }
}


AbstractVar *Match::findVariable(const std::string &name) {
    for(VarList::iterator var = theVariables.begin();
        var != theVariables.end(); ++var) {
        if((*var)->getName() == name) return *var;
    }

    return 0;
}


void Match::getVariables(Vector<double> &x) const {
    x = Vector<double>(theVariables.size(), 0.0);
    int j = 0;

    for(VarList::const_iterator var = theVariables.begin();
        var != theVariables.end(); ++var) {
        x[j++] = (*var)->getInternalValue();
    }
}


void Match::setVariables(const Vector<double> &x) {
    int j = 0;

    for(VarList::iterator var = theVariables.begin();
        var != theVariables.end(); ++var) {
        (*var)->setInternalValue(x[j++]);
    }

    // Force update.
    OpalData::getInstance()->makeDirty(0);
}


int Match::countVariables() const {
    return theVariables.size();
}


void Match::addFunction(AbstractFun *fun) {
    theFunctions.push_back(fun);
    constraintCount += fun->countConstraints();
}


bool Match::evaluate(const Vector<double> &x, Vector<double> &f) {
    // Evaluate the matching functions.
    callCount++;
    setVariables(x);

    // Evaluate and copy the matching function values.
    OpalData::getInstance()->update();
    getFunctions(f);
    return true;
}


void Match::getFunctions(Vector<double> &f) const {
    f = Vector<double>(constraintCount, 0.0);
    int j = 0;

    for(FunList::const_iterator fun = theFunctions.begin();
        fun != theFunctions.end(); ++fun) {
        (*fun)->evaluate(f, j);
    }
}


int Match::countFunctions() const {
    return constraintCount;
}


void Match::print(const char *method, MatchState state) {
    static const char *str[] = {
        "internal",           // Internal call.
        "progress",           // New minimum found.
        "start",              // First call.
        "restart",            // Algorithm restarted.
        "check",              // Check phase.
        "converged",          // Match converged.
        "failed",             // Match failed.
        "call limit",         // Call limit exceeded.
        "accuracy_limit",     // Tolerance too small.
        "terminated"          // Terminated by user.
    };

    // Test if printing desired at all.
    if(( state == INTERNAL && printLevel <= 2 ) ||
       ( state == PROGRESS && printLevel <= 1 )) return;

    // Print variable values.
    cout << "\nMethod: " << method
         << ", Matching state: " << str[state] << ", function calls = "
         << callCount << '\n';

    cout << "\nCurrent variable values:\n";
    for(VarList::const_iterator var = theVariables.begin();
        var != theVariables.end(); ++var) {
        (*var)->print(cout);
    }

    // Print function values.
    Vector<double> f(constraintCount, 0.0);
    callCount++;
    if(( state == INTERNAL && printLevel <= 3 ) ||
       ( state == PROGRESS && printLevel <= 2 )) {
        getFunctions(f);
    } else {
        cout << "\nCurrent matching conditions:\n";
        OpalData::getInstance()->update();
        int index = 0;
        for(FunList::const_iterator fun = theFunctions.begin();
            fun != theFunctions.end(); ++fun) {
            (*fun)->evaluate(f, index);
            (*fun)->print(cout);
        }
    }

    // The sum of squares.
    double sum = 0.0;
    for(Vector<double>::iterator i = f.begin(); i != f.end(); ++i) {
        sum += (*i) * (*i);
    }

    std::streamsize old_prec = cout.precision(12);
    cout << "\nSum of squares = " << sum << std::endl;
    cout.precision(old_prec);
}


int Match::getPrintLevel() const {
    return printLevel;
}


void Match::setPrintLevel(int level) {
    printLevel = level;
}


int Match::getCallCount() const {
    return callCount;
}