// ------------------------------------------------------------------------
// $RCSfile: ConstraintCmd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ConstraintCmd
//   The class for the OPAL CONSTRAINT command.
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/13 15:10:01 $
// $Author: jowett $
//
// ------------------------------------------------------------------------

#include "Match/ConstraintCmd.h"
#include "Attributes/Attributes.h"
#include "AbstractObjects/Expressions.h"
#include "Match/ConcreteFun.h"
#include "Match/Match.h"
#include "Parser/Statement.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/ParseError.h"
#include <iostream>
#include <string>


// Class ConstraintCmd
// ------------------------------------------------------------------------

// The attributes of class ConstraintCmd.
namespace {
    enum {
        LHS,    // the left-hand side of the constraint array
        RHS,    // the right-hand side of the constraint array
        WGT,    // the matching weight array
        SIZE
    };
}


ConstraintCmd::ConstraintCmd():
    Action(SIZE, "CONSTRAINT", "") {
    itsAttr[LHS] = Attributes::makeRealArray
                   ("LHS", "left-hand side");
    itsAttr[RHS] = Attributes::makeRealArray
                   ("RHS", "right-hand side");
    itsAttr[WGT] = Attributes::makeRealArray
                   ("WGT", "the weight(s) for this constraint");
    relation = 0;

    registerOwnership(AttributeHandler::COMMAND);
}


ConstraintCmd::ConstraintCmd(const std::string &name, ConstraintCmd *parent):
    Action(name, parent) {
    relation = 0;
}


ConstraintCmd::~ConstraintCmd()
{}


ConstraintCmd *ConstraintCmd::clone(const std::string &name) {
    return new ConstraintCmd(name, this);
}


void ConstraintCmd::execute() {
    ConcreteFun *fun =
        new ConcreteFun(itsAttr[LHS], relation, itsAttr[RHS], itsAttr[WGT]);
    Match::block->addFunction(fun);

    if(Options::info) {
        std::cerr << "\n 1 constraint has been added.\n" << std::endl;
    }
}


void ConstraintCmd::parse(Statement &stat) {
    Expressions::parseDelimiter(stat, ',');

    // Constraint expressions.
    itsAttr[LHS].parse(stat, false);

    if(stat.delimiter("==")) {
        relation = 0;
    } else if(stat.delimiter('>')) {
        relation = 1;
    } else if(stat.delimiter('<')) {
        relation = 2;
    } else {
        throw ParseError("ConstraintCmd::parse()",
                         "Expected relational operator: "
                         "\"==\", \">\", or \"<\".");
    }

    itsAttr[RHS].parse(stat, false);

    // Weight array.
    Expressions::parseDelimiter(stat, ',');
    std::string name = Expressions::parseString(stat, "Attribute name expected.");

    if(name == "WGT") {
        if(stat.delimiter('=')) {
            itsAttr[WGT].parse(stat, true);
        } else if(stat.delimiter(":=")) {
            throw ParseError("ConstraintCmd::parse()",
                             "The delimiter \":=\" is not allowed here.");
        } else {
            throw ParseError("ConstraintCmd::parse()",
                             "The attribute \"WGT\" has no default.");
        }
    } else {
        throw ParseError("ConstraintCmd::parse()",
                         "Object \"CONSTRAINT\" has no attribute \"" +
                         name + "\".");
    }
}


void ConstraintCmd::print(std::ostream &os) const {
    os << getOpalName();
    Object *parent = getParent();
    if(parent != 0  &&  ! parent->getOpalName().empty()) {
        if(! getOpalName().empty()) os << ':';
        os << parent->getOpalName();
    }

    os << "," << itsAttr[LHS]
       << ((relation == 0) ? "==" : ((relation > 0) ? ">" : "<"))
       << itsAttr[RHS] << ",WGT=" << itsAttr[WGT];
    os << ';';
    os << std::endl;
}


void ConstraintCmd::printHelp(std::ostream &os) const {
    os << "\nThe \"CONSTRAINT\" sub-command defines a constraint.\n"
       << "Its format is:\n"
       << "\tCONSTRAINT,<lhs><relop><rhs>,W=<wgt>"
       << "where <lhs>, <rhs>, and <wgt> are array expressions,"
       << "and <relop> is one of \"<\", \"==\", or \">\"." << std::endl;
}