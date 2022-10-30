#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Classic/AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "Elements/OpalEnge.h"

namespace PyOpal {
namespace PyEnge {

const char* module_docstring = 
    "enge module holds an Enge end field 'plugin' class for use with field models";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalEnge>::attributes = {
    {"X0", "x0", "", PyOpalObjectNS::DOUBLE},
    {"LAMBDA", "enge_lambda", "", PyOpalObjectNS::DOUBLE}, // lambda is a python reserved word
    {"COEFFICIENTS", "coefficients", "", PyOpalObjectNS::FLOAT_LIST},
};

double function(PyOpalObjectNS::PyOpalObject<OpalEnge> pyobject, double x, int n) {
    std::shared_ptr<OpalEnge> objectPtr = pyobject.getOpalShared();
    objectPtr->update();
    std::string name = objectPtr->getOpalName();
    std::shared_ptr<endfieldmodel::EndFieldModel> model = endfieldmodel::EndFieldModel::getEndFieldModel(name);
    model->setMaximumDerivative(n);
    double value = model->function(x, n);
    return value;
}

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalEnge>::classDocstring = 
"Enge class is a field element that models a Enge function.\n"
"\n"
"It is referenced and set to field elements by use of the OPAL name (note not\n"
"any python name). The OPAL name is set using 'set_opal_name'.\n";

BOOST_PYTHON_MODULE(enge) { 
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalEnge> element;
    auto elementClass = element.make_class("Enge");
    elementClass.def("function", &function);
}

}
}
