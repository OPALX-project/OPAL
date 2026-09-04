#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "Elements/OpalPolynomialTimeDependence.h"
#include "PyOpal/PyElements/PyAbstractTimeDependence.h"

namespace PyOpal {
    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalPolynomialTimeDependence>::attributes = {
            {"P0", "p0", "", DOUBLE},
            {"P1", "p1", "", DOUBLE},
            {"P2", "p2", "", DOUBLE},
            {"P3", "p3", "", DOUBLE},
            {"COEFFICIENTS", "coefficients", "", FLOAT_LIST},
    };

    template <>
    std::string PyOpalObjectNS::PyOpalObject<OpalPolynomialTimeDependence>::classDocstring =
        "PolynomialTimeDependence class enables modelling time dependent behaviour.\n"
        "\n"
        "PolynomialTimeDependence class enables modelling lattice parameters that are\n"
        "changing in time, for example RF voltages, frequencies, etc. The dependence\n"
        "is referenced and set to field elements by use of the OPAL name (note not\n"
        "any python name). The OPAL name is set/retrieved using 'set_opal_name' and\n"
        "'get_opal_name' respectively.\n";

    namespace PyOpalPolynomialTimeDependence {

        const char* module_docstring =
            "polynomial_time_dependence module holds a PolynomialTimeDependence class";

        BOOST_PYTHON_MODULE(polynomial_time_dependence) {
            Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalPolynomialTimeDependence> element;
            auto elementClass = element.make_class("PolynomialTimeDependence");
            elementClass.def(
                "function",
                &PyOpal::PyAbstractTimeDependence::function<OpalPolynomialTimeDependence>);
            elementClass.def("update", &PyOpalObjectNS::update<OpalPolynomialTimeDependence>);
        }

    }  // namespace PyOpalPolynomialTimeDependence
}  // namespace PyOpal
