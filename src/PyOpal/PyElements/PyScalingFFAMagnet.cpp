#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Classic/AbsBeamline/ScalingFFAMagnet.h"
#include "Elements/OpalScalingFFAMagnet.h"

namespace PyOpal {

/** Overload this method on the PyOpalObject - it is supposed to be called every
 *  time the magnet is updated.
 */
template <>
void PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet>::doSetup() {
    ScalingFFAMagnet* magnet =
                      dynamic_cast<ScalingFFAMagnet*>(object_m->getElement());
    magnet->setupEndField();
}

namespace PyScalingFFAMagnet {

const char* module_docstring = 
    "scaling_ffa_magnet contains the ScalingFFAMagnet class";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet>::attributes = {
    {"B0", "b0", "", PyOpalObjectNS::DOUBLE},
    {"R0", "r0", "", PyOpalObjectNS::DOUBLE},
    {"FIELD_INDEX", "field_index", "", PyOpalObjectNS::DOUBLE},
    {"TAN_DELTA", "tan_delta", "", PyOpalObjectNS::DOUBLE},
    {"MAX_Y_POWER", "max_vertical_power", "", PyOpalObjectNS::INT},
    {"END_FIELD_MODEL", "end_field_model", "", PyOpalObjectNS::STRING},
    {"END_LENGTH", "end_length", "", PyOpalObjectNS::DOUBLE},
    {"CENTRE_LENGTH", "centre_length", "", PyOpalObjectNS::DOUBLE},
    {"HEIGHT", "height", "", PyOpalObjectNS::DOUBLE},
    {"RADIAL_NEG_EXTENT", "radial_neg_extent", "", PyOpalObjectNS::DOUBLE},
    {"RADIAL_POS_EXTENT", "radial_pos_extent", "", PyOpalObjectNS::DOUBLE},
    {"MAGNET_START", "magnet_start", "", PyOpalObjectNS::DOUBLE},
    {"MAGNET_END", "magnet_end", "", PyOpalObjectNS::DOUBLE},
    {"AZIMUTHAL_EXTENT", "azimuthal_extent", "", PyOpalObjectNS::DOUBLE},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet>::classDocstring = 
"ScalingFFAMagnet class is a field element that models a Scaling FFA magnet.";


BOOST_PYTHON_MODULE(scaling_ffa_magnet) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalScalingFFAMagnet> element;
    auto elementClass = element.make_class("ScalingFFAMagnet");
    element.addGetOpalElement(elementClass);
    element.addGetFieldValue(elementClass, 1e+3, 1.0, 1.0, 1e-1);
}

}
}
