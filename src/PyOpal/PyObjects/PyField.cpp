//
// Python API for Global (external) field access
//
// Copyright (c) 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#include <boost/python.hpp>

#include "Utilities/OpalException.h"
#include "AbsBeamline/Ring.h"
#include "Track/TrackRun.h"
#include "Algorithms/ParallelTTracker.h"
#include "Algorithms/ParallelCyclotronTracker.h"
#include "Physics/Units.h"

#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"

namespace PyOpal {
namespace Field {

std::string field_docstring = 
  "field module enables user to get the field at a point";

std::string get_field_value_docstring =
  "Get the field value at a point in the field map.\n"
  "Only available in Cyclotron mode.\n"
  "\n"
  "The field lookup is performed against the last RINGDEFINITION that was\n"
  "instantiated. This should be instantiated by calling\n"
  "pyopal.parser.initialise_from_opal_file\n"
  "\n"
  "Parameters\n"
  "----------\n"
  "x : float\n"
  "    x position [m]\n"
  "y : float\n"
  "    y position [m]\n"
  "z : float\n"
  "    z position [m]\n"
  "t: float\n"
  "    time [ns]\n"
  "\n"
  "Returns\n"
  "-------\n"
  "The function returns a tuple containing 7 values:\n"
  "out of bounds : int\n"
  "    1 if the event was out of the field map boundary, else 0.\n"
  "Bx : float\n"
  "    x magnetic field [T]\n"
  "By : float\n"
  "    y magnetic field [T]\n"
  "Bz : float\n"
  "    z magnetic field [T]\n"
  "Ex : float\n"
  "    x electric field\n"
  "Ey : float\n"
  "    y electric field\n"
  "Ez : float\n"
  "    z electric field\n";

py::object get_field_value_cyclotron(double x,
                                     double y,
                                     double z,
                                     double t,
                                     ParallelCyclotronTracker* tracker) {
    if (tracker == NULL) {
        throw(OpalException("PyField::get_field_value_cyclotron",
                            "ParallelCyclotronTracker was NULL"));
    }
    Vector_t R(x, y, z);
    Vector_t P, B, E;
    int outOfBounds = tracker->computeExternalFields_m(R, P, t, E, B);
    boost::python::tuple value = boost::python::make_tuple(outOfBounds,
                                          B[0], B[1], B[2], E[0], E[1], E[2]);
    return value;

}

py::object get_field_value(double x, double y, double z, double t) {
    std::shared_ptr<Tracker> tracker = TrackRun::getTracker();
    ParallelCyclotronTracker* trackerCycl = 
                        dynamic_cast<ParallelCyclotronTracker*>(tracker.get());
    if (trackerCycl != nullptr) {
        return get_field_value_cyclotron(x, y, z, t, trackerCycl);
    }
    throw(OpalException("PyField::get_field_value",
                        "Could not find a ParallelCyclotronTracker - get_field_value only works in OPAL-CYCL mode"));
}



// returns a *borrowed* pointer
Ring* getRing() {
    std::shared_ptr<Tracker> tracker = TrackRun::getTracker();
    ParallelCyclotronTracker* trackerCycl = 
                        dynamic_cast<ParallelCyclotronTracker*>(tracker.get());
    Ring* ring = trackerCycl->getRing();
    if (ring == nullptr) {
        throw GeneralClassicException("PyRingDefinition::getSection",
              "Internal PyOpal error - failed to cast to a Ring object");
    } 
    return ring;
}


std::string element_name_docstring = 
"Return a string holding the name of the i^th element [m].\n\n";
std::string getElementName(int i) {
    RingSection* sec = getRing()->getSection(i);
    Component* component = sec->getComponent();
    if (component == nullptr) {
        throw GeneralClassicException("PyRingDefinition::getElementName",
              "Internal PyOpal error - failed to cast to a Component");
    }
    return component->getName();
}

std::string start_pos_docstring = 
"Return a tuple holding the start position of the element (x, y, z) [m].\n\n";

boost::python::object getElementStartPosition(int i) {
    RingSection* sec = getRing()->getSection(i);
    Vector_t pos = sec->getStartPosition();
    return boost::python::make_tuple(pos[0]*Units::mm2m,
                                     pos[1]*Units::mm2m,
                                     pos[2]*Units::mm2m);
}

std::string end_pos_docstring = 
"Return a tuple holding the end position of the element (x, y, z) [m].\n\n";
boost::python::object getElementEndPosition(int i) {
    RingSection* sec = getRing()->getSection(i);
    Vector_t pos = sec->getEndPosition();
    return boost::python::make_tuple(pos[0]*Units::mm2m,
                                     pos[1]*Units::mm2m,
                                     pos[2]*Units::mm2m);
}

std::string start_norm_docstring = 
"Return a tuple holding the vector (x, y, z) normal to the face of the\n"
"element start, pointing towards the element and having length 1.\n\n";
boost::python::object getElementStartNormal(int i) {
    RingSection* sec = getRing()->getSection(i);
    Vector_t dir = sec->getStartNormal();
    return boost::python::make_tuple(dir[0], dir[1], dir[2]);
}

std::string end_norm_docstring = 
"Return a tuple holding the vector (x, y, z) normal to the face of the\n"
"element end, pointing towards the next element and having length 1.\n\n";
boost::python::object getElementEndNormal(int i) {
    RingSection* sec = getRing()->getSection(i);
    Vector_t dir = sec->getEndNormal();
    return boost::python::make_tuple(dir[0], dir[1], dir[2]);
}

std::string num_elements_docstring = 
"Return an integer corresponding to the number of elements stored in the Ring\n"
"If this is 0, check that the track has been executed - the element\n"
"placements are done during Track setup.\n\n";
size_t getNumberOfElements() {
    return getRing()->getNumberOfRingSections();
}


BOOST_PYTHON_MODULE(field) {
    ExceptionTranslation::registerExceptions();
    PyOpal::Globals::Initialise();
    py::scope().attr("__doc__") = field_docstring.c_str();
    py::def("get_field_value",
            get_field_value,
            py::args("x", "y", "z", "t"),
            get_field_value_docstring.c_str()
    );
    py::def("get_number_of_elements",
            getNumberOfElements,
            num_elements_docstring.c_str());
    py::def("get_element_start_position",
            getElementStartPosition,
            element_name_docstring.c_str());
    py::def("get_element_start_normal", &getElementStartNormal);
    py::def("get_element_end_position", &getElementEndPosition);
    py::def("get_element_end_normal", &getElementEndNormal);
    py::def("get_element_name", &getElementName);
}

}
}

