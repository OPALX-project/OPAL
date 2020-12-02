#include <Python.h>
#include <structmember.h>

#include <exception>
#include <iostream>
#include <boost/python.hpp>

#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "AbsBeamline/Ring.h"
#include "Elements/OpalRingDefinition.h"

//using namespace boost::python;

std::string get_field_value_docstring = 
std::string("Get the field value at a point in the field map.\n\n")+
std::string("    x: x position [m]\n")+
std::string("    y: y position [m]\n")+
std::string("    z: z position [m]\n")+
std::string("    t: time [ns]\n")+
std::string("Returns a tuple containing 6 values:\n")+
std::string("    out of bounds: 1 if the event was out of the field map\n")+
std::string("                   boundary, else 0.\n")+
std::string("    Bx: x magnetic field [T]\n")+
std::string("    By: y magnetic field [T]\n")+
std::string("    Bz: z magnetic field [T]\n")+
std::string("    Ex: x electric field\n")+
std::string("    Ey: y electric field\n")+
std::string("    Ez: z electric field\n");

boost::python::object get_field_value(double x, double y, double z, double t) {
    Ring* ring = const_cast<Ring*>(Ring::getLastLockedRing());
    if (ring == NULL) {
        std::string err = "Could not find a ring object - maybe a "+
           std::string("RingDefinition was not defined or KeepAlive was False");
        // need to define proper exception handling
        throw(err);
    }
    Vector_t R(x, y, z);
    Vector_t P(0, 0, 0);
    Vector_t E, B;
    int outOfBounds = ring->apply(R, P, t, E, B);
    boost::python::tuple value = boost::python::make_tuple(outOfBounds,
                                          B[0]/10., B[1]/10., B[2]/10.,
                                          E[0], E[1], E[2]);
    return value;
}

const char* module_docstring = "field module returns the field";

BOOST_PYTHON_MODULE(boost_field) {
    boost::python::def("get_field_value",
                       get_field_value,
                       boost::python::args("x", "y", "z", "t"),
                       get_field_value_docstring.c_str());
}
