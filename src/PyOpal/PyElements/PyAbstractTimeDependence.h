#ifndef PYELEMENTS_PYABSTRACTTIMEDEPENDENCE_H
#define PYELEMENTS_PYABSTRACTTIMEDEPENDENCE_H

#include "PyOpal/PyCore/PyOpalObject.h"
#include "Classic/Algorithms/AbstractTimeDependence.h"

namespace PyOpal {
namespace PyAbstractTimeDependence {

std::string functionDocstring = 
"Update the TimeDependence and then calculate the value at a give time\n"
"- t [float]: time at which to evaluate the function\n"
"Returns the function value.\n";

template <class C>
double function(PyOpalObjectNS::PyOpalObject<C> pyobject, double t) {
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    objectPtr->update();
    std::string name = objectPtr->getOpalName();
    std::shared_ptr<AbstractTimeDependence> model =
                                AbstractTimeDependence::getTimeDependence(name);
    double value = model->getValue(t);
    return value;
}

}
}

#endif // PYELEMENTS_PYABSTRACTTIMEDEPENDENCE_H
