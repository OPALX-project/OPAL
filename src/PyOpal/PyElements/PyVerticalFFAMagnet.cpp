#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Elements/OpalVerticalFFAMagnet.h"

namespace PyOpal {
namespace PyVerticalFFAMagnet {

const char* module_docstring = 
    "vertical_ffa_magnet contains the VerticalFFAMagnet class";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet>::attributes = {
    {"B0", "b0", "", PyOpalObjectNS::DOUBLE},
    {"FIELD_INDEX", "field_index", "", PyOpalObjectNS::DOUBLE},
    {"WIDTH", "width", "", PyOpalObjectNS::DOUBLE},
    {"MAX_HORIZONTAL_POWER", "max_horizontal_power", "", PyOpalObjectNS::INT},
    {"END_LENGTH", "end_length", "", PyOpalObjectNS::DOUBLE},
    {"CENTRE_LENGTH", "centre_length", "", PyOpalObjectNS::DOUBLE},
    {"BB_LENGTH", "bb_length", "", PyOpalObjectNS::DOUBLE},
    {"HEIGHT_NEG_EXTENT", "height_neg_extent", "", PyOpalObjectNS::DOUBLE},
    {"HEIGHT_POS_EXTENT", "height_pos_extent", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet>::classDocstring = 
"VerticalFFAMagnet class is a field element that models a Vertical FFA magnet.";

BOOST_PYTHON_MODULE(vertical_ffa_magnet) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalVerticalFFAMagnet> element;
    auto elementClass = element.make_class("VerticalFFAMagnet");
    element.addGetFieldValue(elementClass, 1e+3, 1.0, 1.0, 1e-1);
}

}
}
