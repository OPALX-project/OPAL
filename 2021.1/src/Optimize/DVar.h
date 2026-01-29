#ifndef OPAL_DVAR_HH
#define OPAL_DVAR_HH

#include "AbstractObjects/Definition.h"

class DVar: public Definition {
public:
    DVar();
    ~DVar();

    virtual DVar *clone(const std::string &name);

    virtual void execute();

    std::string getVariable() const;
    double getLowerBound() const;
    double getUpperBound() const;

private:
    DVar(const std::string &name,
         DVar *parent);
};

inline
DVar* DVar::clone(const std::string &name) {
    return new DVar(name, this);
}

#endif