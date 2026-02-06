#include "AbsBeamline/EndFieldModel/AsymmetricEnge.h"

namespace endfieldmodel {

AsymmetricEnge::AsymmetricEnge() : engeStart_m(new Enge()),
    engeEnd_m(new Enge()) {
}

AsymmetricEnge::AsymmetricEnge(const AsymmetricEnge& rhs) 
  : engeStart_m(rhs.engeStart_m->clone()), engeEnd_m(rhs.engeEnd_m->clone()) {
}

AsymmetricEnge::AsymmetricEnge(const std::vector<double> aStart,
                       double x0Start,
                       double lambdaStart, 
                       const std::vector<double> aEnd,
                       double x0End,
                       double lambdaEnd) : engeStart_m(new Enge()),
    engeEnd_m(new Enge()) {
    engeStart_m->setCoefficients(aStart);
    engeStart_m->setX0(x0Start);
    engeStart_m->setLambda(lambdaStart);
    // x0 is held in this
    engeEnd_m->setCoefficients(aEnd);
    engeEnd_m->setX0(x0End);
    engeEnd_m->setLambda(lambdaEnd);
}

void AsymmetricEnge::rescale(double scaleFactor) {
    engeStart_m->rescale(scaleFactor);
    engeEnd_m->rescale(scaleFactor);
}

std::ostream& AsymmetricEnge::print(std::ostream& out) const {
    out << "AsymmetricEnge start ";
    engeStart_m->print(out);
    out  << " end ";
    engeStart_m->print(out);
    return out;
}


double AsymmetricEnge::function(double x, int n) const {
    // f(x) = E(x-x0) + E(-x-x0) - 1
    // f^{(2n)} = E^{(2n)}(x-x0) + E^{(2n)}(-x-x0)
    // f^{(2n+1)} = E^{(2n)}(x-x0) - E^{(2n)}(-x-x0)
    double sup = -x-engeStart_m->getX0();
    double sdn = +x-engeEnd_m->getX0();
    if (n == 0) {
        return engeStart_m->getEnge(sup, n)+
               engeEnd_m->getEnge(sdn, n)-1;
    } else if (n%2) {
        return -engeStart_m->getEnge(sup, n)+
               engeEnd_m->getEnge(sdn, n);
    } else {
        return engeStart_m->getEnge(sup, n)+
               engeEnd_m->getEnge(sdn, n);
    }
}

}
