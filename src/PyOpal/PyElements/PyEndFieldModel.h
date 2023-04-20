#ifndef PYELEMENTS_PYENDFIELDMODEL_H
#define PYELEMENTS_PYENDFIELDMODEL_H

#include "PyOpal/PyCore/PyOpalObject.h"
#include "Classic/AbsBeamline/EndFieldModel/EndFieldModel.h"

namespace PyOpal {
namespace PyEndFieldModel {

std::string functionDocstring = 
"update the Enge object and then calculate the enge function or derivative\n"
"- x [float]: position at which to evaluate the function\n"
"- n [int]: derivative to calculate; 0 will return the enge function, 1 will\n"
"           return the 1st derivative, etc. Note that floating point precision\n"
"           becomes increasingly limiting for n > 10 or so.\n"
"Returns the function value or derivative.\n";

template <class C>
double function(PyOpalObjectNS::PyOpalObject<C> pyobject, double x, int n) {
    if (n < 0) {
        throw OpalException("n must be positive as it indexes the derivative", 
                            "PyEndFieldModel::function");
    }
    std::shared_ptr<C> objectPtr = pyobject.getOpalShared();
    objectPtr->update();
    std::string name = objectPtr->getOpalName();
    std::shared_ptr<endfieldmodel::EndFieldModel> model = endfieldmodel::EndFieldModel::getEndFieldModel(name);
    model->setMaximumDerivative(n);
    double value = model->function(x, n);
    return value;
}

}
}

#endif // PYELEMENTS_PYENDFIELDMODEL_H
