#include <Python.h>
#include <structmember.h>
#include <iostream>

#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/Object.h"
#include "AbsBeamline/Ring.h"
#include "Elements/OpalRingDefinition.h"


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

PyObject* get_field_value(PyObject */*self*/, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = {
        const_cast<char*>("x"),
        const_cast<char*>("y"),
        const_cast<char*>("z"),
        const_cast<char*>("t"),
        NULL};
    double x, y, z, t;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
                                    &x, &y, &z, &t)) {
        return NULL;
    }
    Ring* ring = const_cast<Ring*>(Ring::getLastLockedRing());
    if (ring == NULL) {
        std::string err = "Could not find a ring object - maybe a "+
           std::string("RingDefinition was not defined or KeepAlive was False");
        PyErr_SetString(PyExc_ValueError, err.c_str());
        return 0;
    }
    Vector_t R(x, y, z);
    Vector_t P(0, 0, 0);
    Vector_t E, B;
    int outOfBounds = ring->apply(R, P, t, E, B);
    PyObject* value = Py_BuildValue("idddddd", outOfBounds,
                                    B[0]/10., B[1]/10., B[2]/10.,
                                    E[0], E[1], E[2]);
    return value;
}


std::string get_field_derivative_docstring = 
std::string("Get the field derivative at a point in the field map.\n\n")+
std::string("    x: x position [m]\n")+
std::string("    y: y position [m]\n")+
std::string("    z: z position [m]\n")+
std::string("    t: time [ns]\n")+
std::string("    field_index: index controlling the field element whose\n")+
std::string("       derivative will be calculated, as element of\n")+
std::string("       dF = (dBx, dBy, dBz, dEx, dEy, dEz)\n")+
std::string("    pos_index: index controlling the position element as an\n")+
std::string("       element of du = (dx, dy, dz, dt)")+
std::string("    delta: position scale [m or ns]\n")+
std::string("Returns the derivative dF/du. The derivative is calculated by\n")+
std::string("dF/du = ((u+du)-(u-du))[pos_index]/\n")+
std::string("        (F(u+du)-F(u-du))[field_index]\n");

PyObject* get_field_derivative(PyObject */*self*/, PyObject *args, PyObject *kwds) {
    std::string err = "Returning fishy results - disabled";
    PyErr_SetString(PyExc_ValueError, err.c_str());
    return 0;
    static char *kwlist[] = {
        const_cast<char*>("x"),
        const_cast<char*>("y"),
        const_cast<char*>("z"),
        const_cast<char*>("t"),
        const_cast<char*>("field_index"),
        const_cast<char*>("pos_index"),
        const_cast<char*>("delta"),
        NULL};
    double x, y, z, t, delta;
    int field_index, pos_index;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ddddiid", kwlist,
                        &x, &y, &z, &t, &field_index, &pos_index, &delta)) {
        return NULL;
    }
    Ring* ring = const_cast<Ring*>(Ring::getLastLockedRing());
    if (ring == NULL) {
        std::string err = "Could not find a ring object - maybe a "+
           std::string("RingDefinition was not defined or KeepAlive was False");
        PyErr_SetString(PyExc_ValueError, err.c_str());
        return 0;
    }
    if (pos_index < 0 or pos_index > 3) {
        std::string err = "pos_index was out of range";
        PyErr_SetString(PyExc_ValueError, err.c_str());
        return 0;
    }
    if (field_index < 0 or field_index > 5) {
        std::string err = "field_index was out of range";
        PyErr_SetString(PyExc_ValueError, err.c_str());
        return 0;
    }
    // now calculate the derivative
    Vector_t deltaR(0., 0., 0.);
    double deltaT = 0;
    if (pos_index < 3) {
        deltaR[pos_index] = delta;
    } else if (pos_index == 3) {
        deltaT = delta;
    }
    Vector_t R = Vector_t(x, y, z)+deltaR;
    Vector_t P(0, 0, 0);
    t += deltaT;
    Vector_t Eupper, Bupper;
    ring->apply(R, P, t, Eupper, Bupper);
    R = R - 2*deltaR;
    t -= 2*deltaT;
    Vector_t Elower, Blower;
    ring->apply(R, P, t, Eupper, Bupper);
    double deltaF = 0.0;
    if (field_index < 3) {
        deltaF = Bupper[field_index]-Blower[field_index];
    } else {
        deltaF = Eupper[field_index-3]-Elower[field_index-3];
    }
    double derivative = deltaF/2/delta;
    PyObject* pyDerivative = PyFloat_FromDouble(derivative);
    Py_INCREF(pyDerivative);
    return pyDerivative;
}


std::string get_potential_docstring = 
std::string("Get the potential at a point in the field map.\n\n")+
std::string("    x: x position [m]\n")+
std::string("    y: y position [m]\n")+
std::string("    z: z position [m]\n")+
std::string("    t: time [ns]\n")+
std::string("Returns a tuple containing 4 values:\n")+
std::string("    out of bounds: 1 if the event was out of the field map\n")+
std::string("                   boundary, else 0.\n")+
std::string("    Ax: x vector potential [Tm]\n")+
std::string("    Ay: y vector potential [Tm]\n")+
std::string("    Az: z vector potential [Tm]\n")+
std::string("    phi: electric potential\n");

PyObject* get_potential(PyObject */*self*/, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = {
        const_cast<char*>("x"),
        const_cast<char*>("y"),
        const_cast<char*>("z"),
        const_cast<char*>("t"),
        NULL};
    double x, y, z, t;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dddd", kwlist,
                                    &x, &y, &z, &t)) {
        return NULL;
    }
    Ring* ring = const_cast<Ring*>(Ring::getLastLockedRing());
    if (ring == NULL) {
        std::string err = "Could not find a ring object - maybe PyOpal.parser "+
           std::string("was not called, a RingDefinition was not set up ")+
           std::string("or KeepAlive was False");
        PyErr_SetString(PyExc_ValueError, err.c_str());
        return 0;
    }
    Vector_t R(x, y, z);
    Vector_t A;
    double phi;
    int outOfBounds = ring->getPotential(R, t, A, phi);
    PyObject* value = Py_BuildValue("idddd", outOfBounds,
                                    A[0], A[1], A[2], phi);
    return value;
}


const char* module_docstring = "field module returns the field";

static PyMethodDef _module_methods[] = {
{"get_field_value", (PyCFunction)get_field_value,
  METH_VARARGS|METH_KEYWORDS, get_field_value_docstring.c_str()},
{"get_field_derivative", (PyCFunction)get_field_derivative,
  METH_VARARGS|METH_KEYWORDS, get_field_value_docstring.c_str()},
{"get_potential", (PyCFunction)get_potential,
  METH_VARARGS|METH_KEYWORDS, get_potential_docstring.c_str()},
{NULL, NULL, 0, NULL}
};


static struct PyModuleDef fielddef = {
    PyModuleDef_HEAD_INIT,
    "field",     /* m_name */
    module_docstring,  /* m_doc */
    -1,                  /* m_size */
    _module_methods,    /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_field(void) {
    //PyOpal::Globals::Initialise();
    PyObject* module = PyModule_Create(&fielddef);
    return module;
}
