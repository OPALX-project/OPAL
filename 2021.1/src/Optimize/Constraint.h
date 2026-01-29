#ifndef OPAL_CONSTRAINT_HH
#define OPAL_CONSTRAINT_HH

#include "AbstractObjects/Definition.h"

class Constraint: public Definition {
public:
    Constraint();
    ~Constraint();

    virtual Constraint *clone(const std::string &name);

    virtual void execute();

    std::string getExpression() const;
private:
    Constraint(const std::string &name,
         Constraint *parent);
};

inline
Constraint* Constraint::clone(const std::string &name) {
    return new Constraint(name, this);
}

#endif