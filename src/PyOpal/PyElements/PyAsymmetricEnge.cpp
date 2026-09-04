#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyElements/PyEndFieldModel.h"

#include "Classic/AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "Elements/OpalAsymmetricEnge.h"

namespace PyOpal {
    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalAsymmetricEnge>::attributes = {
            {"X0_START", "x0_start", "", PyOpalObjectNS::DOUBLE},
            {"LAMBDA_START", "lambda_start", "", PyOpalObjectNS::DOUBLE},
            {"COEFFICIENTS_START", "coefficients_start", "", PyOpalObjectNS::FLOAT_LIST},
            {"X0_END", "x0_end", "", PyOpalObjectNS::DOUBLE},
            {"LAMBDA_END", "lambda_end", "", PyOpalObjectNS::DOUBLE},
            {"COEFFICIENTS_END", "coefficients_end", "", PyOpalObjectNS::FLOAT_LIST},
    };

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalAsymmetricEnge>::classDocstring =
        "AsymmetricEnge class is a field plugin that models a AsymmetricEnge function.\n"
        "\n"
        "It is referenced and set to field elements by use of the OPAL name (note not\n"
        "any python name). The OPAL name is set using 'set_opal_name'.\n";

    namespace PyAsymmetricEnge {

        const char* module_docstring =
            "asymmetric_enge module holds an AsymmetricEnge end field 'plugin' class for\n"
            "use with field models\n";

        BOOST_PYTHON_MODULE(asymmetric_enge) {
            PyOpal::Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalAsymmetricEnge> element;
            auto elementClass = element.make_element_class("AsymmetricEnge");
            elementClass.def("function", &PyOpal::PyEndFieldModel::function<OpalAsymmetricEnge>);
            elementClass.def("update", &PyOpalObjectNS::update<OpalAsymmetricEnge>);
        }

    }  // namespace PyAsymmetricEnge
}  // namespace PyOpal
