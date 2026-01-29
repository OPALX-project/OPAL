#include "Optimize/Objective.h"
#include "Attributes/Attributes.h"

namespace {
    enum {
        EXPR,
        SIZE
    };
}

Objective::Objective():
    Definition(SIZE, "OBJECTIVE", "The OBJECTIVE statement defines an objective  for optimization")
{
    itsAttr[EXPR] = Attributes::makeString("EXPR",
                                           "Expression to minimize during optimization");

    registerOwnership(AttributeHandler::STATEMENT);
}

Objective::Objective(const std::string &name, Objective *parent):
    Definition(name, parent)
{ }

Objective::~Objective()
{ }

void Objective::execute() {

}

std::string Objective::getExpression() const {
    std::string expr = Attributes::getString(itsAttr[EXPR]);

    return expr;
}