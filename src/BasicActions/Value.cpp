// ------------------------------------------------------------------------
// $RCSfile: Value.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Value
//   The class for the OPAL VALUE command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:37 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "BasicActions/Value.h"

#include "Attributes/Attributes.h"
#include "Parser/Statement.h"
#include "Utility/Inform.h"

#include <string>
#include <vector>


// Class Value
// ------------------------------------------------------------------------

extern Inform *gmsg;

Value::Value():
    Action(1, "VALUE",
           "The \"VALUE\" statement prints a list of expressions and "
           "their values.") {
    itsAttr[0] = Attributes::makeRealArray
                 ("VALUE", "The values to be evaluated");

    registerOwnership(AttributeHandler::STATEMENT);
}


Value::Value(const std::string &name, Value *parent):
    Action(name, parent)
{}


Value::~Value()
{}


Value *Value::clone(const std::string &name) {
    return new Value(name, this);
}


void Value::execute() {
    *gmsg << "\nvalue: " << itsAttr[0] << " = {";
    //  std::streamsize old_prec = *gmsg.precision(12);
    const std::vector<double> array = Attributes::getRealArray(itsAttr[0]);
    std::vector<double>::const_iterator i = array.begin();

    while(i != array.end()) {
        *gmsg << *i++;
        if(i == array.end()) break;
        *gmsg << ", ";
    }

    *gmsg << "}\n" << endl;
    //  *gmsg.precision(old_prec);
}

void Value::parse(Statement &statement) {
    // parse, but don't evaluate (for printing mainly)
    parseShortcut(statement, false);
}