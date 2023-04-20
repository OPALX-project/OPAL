#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Elements/OpalOffset/OpalLocalCartesianOffset.h"

//using namespace boost::python;
namespace PyOpal {
namespace PyOpalLocalCartesianOffset {

using OpalOffset::OpalLocalCartesianOffset;

const char* module_docstring = "build a local cartesian offset";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalLocalCartesianOffset>::attributes = {
    {"END_POSITION_X", "end_position_x", "", PyOpalObjectNS::DOUBLE},
    {"END_POSITION_Y", "end_position_y", "", PyOpalObjectNS::DOUBLE},
    {"END_NORMAL_X", "end_normal_x", "", PyOpalObjectNS::DOUBLE},
    {"END_NORMAL_Y", "end_normal_y", "", PyOpalObjectNS::DOUBLE},
};

BOOST_PYTHON_MODULE(local_cartesian_offset) {
    Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalLocalCartesianOffset> element;
    auto elementClass = element.make_element_class("LocalCartesianOffset");
}

}
}
