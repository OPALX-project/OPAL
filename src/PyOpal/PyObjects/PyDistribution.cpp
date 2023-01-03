#include "Distribution/Distribution.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"

namespace PyOpal {
namespace PyDistributionNS {

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<Distribution>::attributes = {
    {"TYPE", "type", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"FNAME", "filename", "", PyOpalObjectNS::STRING},
    {"INPUTMOUNITS", "momentum_units", "", PyOpalObjectNS::PREDEFINED_STRING},
};

void registerDistribution(PyOpalObjectNS::PyOpalObject<Distribution>& dist) {
    Object* obj = &(*dist.getOpalShared());
    obj->update();
    OpalData::getInstance()->define(obj);
}

BOOST_PYTHON_MODULE(distribution) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<Distribution> distributionObject;
    auto distributionClass = distributionObject.make_class("Distribution");
    distributionObject.addExecute(distributionClass);
    distributionClass.def("register", &registerDistribution);

}

} // PyDistribution
} // PyOpal

