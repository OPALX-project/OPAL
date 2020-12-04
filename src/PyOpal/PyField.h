#include <boost/python.hpp>

#include "PyOpal/Globals.h"

namespace PyOpal {
namespace Field {
std::string field_docstring = 
  "field module enables user to get the field at a point";

std::string get_field_value_docstring =
  "Get the field value at a point in the field map.\n"
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
  "The function returns a tuple containing 6 values:\n"
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

py::object get_field_value(double x, double y, double z, double t);

}
}