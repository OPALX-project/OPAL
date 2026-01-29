#include "Optimize/Constraint.h"
#include "Attributes/Attributes.h"

namespace {
    enum {
        EXPR,
        SIZE
    };
}

Constraint::Constraint():
    Definition(SIZE, "CONSTRAINT", "The CONSTRAINT statement defines a constraint for optimization")
{
    itsAttr[EXPR] = Attributes::makeString("EXPR",
                                           "Expression that should be fulfilled during optimization");

    registerOwnership(AttributeHandler::STATEMENT);
}

Constraint::Constraint(const std::string &name, Constraint *parent):
    Definition(name, parent)
{ }

Constraint::~Constraint()
{ }

void Constraint::execute() {

}

std::string Constraint::getExpression() const {
    std::string expr = Attributes::getString(itsAttr[EXPR]);

    return expr;
}