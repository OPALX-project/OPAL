#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/PyOpalObject.h"

#include "AbstractObjects/OpalData.h"
#include "Track/TrackRun.h"


extern Inform *gmsg;

namespace PyOpal {
namespace PyTrackRunNS {

std::string track_run_docstring = std::string();

const char* module_docstring = "build a tracking object";

template <>
std::vector<PyOpalObjectNS::AttributeDef> PyOpalObjectNS::PyOpalObject<TrackRun>::attributes = {
    {"METHOD", "method", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"TURNS", "turns", "", PyOpalObjectNS::DOUBLE},
    {"MBMODE", "multibunch_mode", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"PARAMB", "multibunch_control", "", PyOpalObjectNS::DOUBLE},
    {"MB_ETA", "multibunch_scale", "", PyOpalObjectNS::DOUBLE},
    {"MB_BINNING", "multibunch_binning", "", PyOpalObjectNS::PREDEFINED_STRING},
    {"BEAM", "beam_name", "", PyOpalObjectNS::STRING},
    {"FIELDSOLVER", "field_solver", "", PyOpalObjectNS::STRING},
    {"BOUNDARYGEOMETRY", "boundary_geometry", "", PyOpalObjectNS::STRING},
    {"DISTRIBUTION", "distribution", "", PyOpalObjectNS::STRING_LIST},
};

template <>
std::string PyOpalObjectNS::PyOpalObject<TrackRun>::classDocstring = "";

void setRunName(PyOpalObjectNS::PyOpalObject<TrackRun>& /*run*/, std::string name) {

    OpalData::getInstance()->storeInputFn(name);
}

BOOST_PYTHON_MODULE(track_run) {
    ExceptionTranslation::registerExceptions();
    PyOpal::Globals::Initialise();
    PyOpalObjectNS::PyOpalObject<TrackRun> trackRun;
    auto trackClass = trackRun.make_class("TrackRun");
    trackRun.addExecute(trackClass);
    trackClass.def("set_run_name", &setRunName);
    setRunName(trackRun, "PyOpal"); // force default run name to "PyOpal"
}

} // PyTrackRun
} // PyOpal

