// ------------------------------------------------------------------------
// $RCSfile: Macro.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Macro/
//   Abstract base class for OPAL macro-like commands.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:43 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "OpalParser/Macro.h"
#include "AbstractObjects/Expressions.h"
#include "Parser/Statement.h"
#include "Utilities/ParseError.h"


// Class Macro
// ------------------------------------------------------------------------

Macro::Macro(int size, const char *name, const char *help):
    Object(size, name, help), formals(), actuals()
{}


Macro::Macro(const std::string &name, Object *parent):
    Object(name, parent), formals(), actuals()
{}


Macro::~Macro()
{}


Macro *Macro::clone(const std::string &) {
    throw ParseError("Macro::clone()",
                     "You cannot execute the command \"" + getOpalName() +
                     "\" without parameters.");
}


const std::string Macro::getCategory() const {
    return "MACRO";
}


bool Macro::shouldTrace() const {
    return false;
}


bool Macro::shouldUpdate() const {
    return false;
}


void Macro::parseActuals(Statement &stat) {
    // We start after the opening '('.
    actuals.clear();

    if(! stat.delimiter(')')) {
        int level = 1;
        std::vector<Token> act;

        while(true) {
            Token token = stat.getCurrent();

            if(token.isDel(',')) {
                // Comma at level=1 terminates an actual argument.
                if(level == 1) {
                    actuals.push_back(act);
                    act.clear();
                    continue;
                }
            } else if(token.isDel('(')) {
                level++;
            } else if(token.isDel(')')) {
                // Closing parenthesis at level=1 terminates the actuals list.
                if(--level == 0) {
                    actuals.push_back(act);
                    act.clear();
                    break;
                }
            }

            // In all other casses append token to current actual.
            act.push_back(token);
        }
    }
}


void Macro::parseFormals(Statement &stat) {
    // We start after the opening '('.
    formals.clear();
    if(! stat.delimiter(')')) {
        do {
            std::string form =
                Expressions::parseString(stat, "Expected formal argument name.");
            formals.push_back(form);
        } while(stat.delimiter(','));

        Expressions::parseDelimiter(stat, ')');
    }

    Expressions::parseDelimiter(stat, ':');
}
