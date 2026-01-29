#include "Optimize/DVar.h"
#include "Attributes/Attributes.h"

namespace {
    enum {
        VARIABLE,
        LOWERBOUND,
        UPPERBOUND,
        SIZE
    };
}

DVar::DVar():
    Definition(SIZE, "DVAR", "The DVAR statement defines a variable for optimization")
{
    itsAttr[VARIABLE] = Attributes::makeString("VARIABLE",
                                               "Variable name that should be varied during optimization");
    itsAttr[LOWERBOUND] = Attributes::makeReal("LOWERBOUND",
                                               "Lower limit of the range of values that the variable should assume");
    itsAttr[UPPERBOUND] = Attributes::makeReal("UPPERBOUND",
                                               "Upper limit of the range of values that the variable should assume");

    registerOwnership(AttributeHandler::STATEMENT);
}

DVar::DVar(const std::string &name, DVar *parent):
    Definition(name, parent)
{ }

DVar::~DVar()
{ }

void DVar::execute() {

}

std::string DVar::getVariable() const {
    return Attributes::getString(itsAttr[VARIABLE]);
}

double DVar::getLowerBound() const {
    return Attributes::getReal(itsAttr[LOWERBOUND]);
}

double DVar::getUpperBound() const {
    return Attributes::getReal(itsAttr[UPPERBOUND]);
}