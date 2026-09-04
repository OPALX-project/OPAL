#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "Elements/OpalOffset/OpalGlobalCartesianOffset.h"

// using namespace boost::python;
namespace PyOpal {
    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<OpalOffset::OpalGlobalCartesianOffset>::attributes = {
            {"END_POSITION_X", "end_position_x", "", PyOpalObjectNS::DOUBLE},
            {"END_POSITION_Y", "end_position_y", "", PyOpalObjectNS::DOUBLE},
            {"END_NORMAL_X", "end_normal_x", "", PyOpalObjectNS::DOUBLE},
            {"END_NORMAL_Y", "end_normal_y", "", PyOpalObjectNS::DOUBLE},
    };

    namespace PyOpalGlobalCartesianOffset {

        using OpalOffset::OpalGlobalCartesianOffset;

        const char* module_docstring = "build a local cartesian offset";

        BOOST_PYTHON_MODULE(global_cartesian_offset) {
            Globals::Initialise();
            ExceptionTranslation::registerExceptions();
            PyOpalObjectNS::PyOpalObject<OpalGlobalCartesianOffset> element;
            auto elementClass = element.make_element_class("GlobalCartesianOffset");
        }

    }  // namespace PyOpalGlobalCartesianOffset
}  // namespace PyOpal
