#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Classic/AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "Elements/OpalAsymmetricEnge.h"

namespace PyOpal {
namespace PyAsymmetricEnge {

const char* module_docstring = 
    "enge module holds an Enge end field 'plugin' class for use with field models";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalAsymmetricEnge>::attributes = {
    {"X0_START", "x0_start", "", PyOpalObjectNS::DOUBLE},
    {"LAMBDA_START", "lambda_start", "", PyOpalObjectNS::DOUBLE},
    {"COEFFICIENTS_START", "coefficients_start", "", PyOpalObjectNS::FLOAT_LIST},
    {"X0_END", "x0_end", "", PyOpalObjectNS::DOUBLE},
    {"LAMBDA_END", "lambda_end", "", PyOpalObjectNS::DOUBLE},
    {"COEFFICIENTS_END", "coefficients_end", "", PyOpalObjectNS::FLOAT_LIST},
};

/** Note very similar code in PyEnge; would be better to share if it is needed much more */
double function(PyOpalObjectNS::PyOpalObject<OpalAsymmetricEnge> pyobject, double x, int n) {
    std::shared_ptr<OpalAsymmetricEnge> objectPtr = pyobject.getOpalShared();
    objectPtr->update();
    std::string name = objectPtr->getOpalName();
    std::shared_ptr<endfieldmodel::EndFieldModel> model = endfieldmodel::EndFieldModel::getEndFieldModel(name);
    model->setMaximumDerivative(n);
    double value = model->function(x, n);
    return value;
}

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalAsymmetricEnge>::classDocstring = 
"AsymmetricEnge class is a field plugin that models a AsymmetricEnge function.\n"
"\n"
"It is referenced and set to field elements by use of the OPAL name (note not\n"
"any python name). The OPAL name is set using 'set_opal_name'.\n";

BOOST_PYTHON_MODULE(asymmetric_enge) { 
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalAsymmetricEnge> element;
    auto elementClass = element.make_class("AsymmetricEnge");
    elementClass.def("function", &function);
}

}
}
