#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/PyOpalObject.h"
#include "PyOpal/PyCore/Globals.h"

#include "Classic/AbsBeamline/Probe.h"
#include "Elements/OpalProbe.h"

namespace PyOpal {
namespace PyProbe {

const char* module_docstring = 
    "probe contains the Probe class";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<OpalProbe>::attributes = {
    {"XSTART", "x_start", "", PyOpalObjectNS::DOUBLE},
    {"XEND", "x_end", "", PyOpalObjectNS::DOUBLE},
    {"YSTART", "y_start", "", PyOpalObjectNS::DOUBLE},
    {"YEND", "y_end", "", PyOpalObjectNS::DOUBLE},
    {"WIDTH", "width", "", PyOpalObjectNS::DOUBLE},
    {"STEP", "step", "", PyOpalObjectNS::DOUBLE},
    {"OUTFN", "output_filename", "", PyOpalObjectNS::STRING}, // OUTFN comes from OpalElement (yes, all elements can have a filename!)
};

template <>
std::string PyOpalObjectNS::PyOpalObject<OpalProbe>::classDocstring = 
"Probe is used to generate output data based on particle tracks crossing a plane.";

BOOST_PYTHON_MODULE(probe) {
    PyOpal::Globals::Initialise();
    ExceptionTranslation::registerExceptions();
    PyOpalObjectNS::PyOpalObject<OpalProbe> element;
    auto elementClass = element.make_element_class("Probe");
}

}
}
