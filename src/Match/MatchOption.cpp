// ------------------------------------------------------------------------
// $RCSfile: MatchOption.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MatchOption
//   The class for the OPAL OPTION command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Match/MatchOption.h"
#include "Attributes/Attributes.h"
#include "Match/Match.h"
#include <cmath>

// Class MatchOption
// ------------------------------------------------------------------------

MatchOption::MatchOption():
    Action(1, "OPTION",
           "The \"OPTION\" sub-command sets matching options.") {
    itsAttr[0] = Attributes::makeReal
                 ("LEVEL", "The desired verbosity for output", 0.0);

    registerOwnership(AttributeHandler::SUB_COMMAND);
}


MatchOption::MatchOption(const std::string &name, MatchOption *parent):
    Action(name, parent)
{}


MatchOption::~MatchOption()
{}


MatchOption *MatchOption::clone(const std::string &name) {
    return new MatchOption(name, this);
}


void MatchOption::execute() {
    int level = int(std::round(Attributes::getReal(itsAttr[0])));
    Match::block->setPrintLevel(level);
}