#ifndef OPAL_OBJECTIVE_HH
#define OPAL_OBJECTIVE_HH

#include "AbstractObjects/Definition.h"

class Objective: public Definition {
public:
    Objective();
    ~Objective();

    virtual Objective *clone(const std::string &name);

    virtual void execute();

    std::string getExpression() const;

private:
    Objective(const std::string &name,
         Objective *parent);
};

inline
Objective* Objective::clone(const std::string &name) {
    return new Objective(name, this);
}

#endif