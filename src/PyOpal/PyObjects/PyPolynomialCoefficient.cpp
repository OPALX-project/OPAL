//
// Python API for PolynomialCoefficient (part of the multidimensional polynomial fitting routines)
//
// Copyright (c) 2008-2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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

#include <Python.h>
#include <structmember.h>

#include "Fields/Interpolation/PolynomialCoefficient.h"
#include "PyOpal/Globals.h"
#include "PyOpal/PyPolynomialCoefficient.h"

namespace PyPolynomialCoefficient {
/*
    PolynomialCoefficient(std::vector<int> inVariablesByVector,
                          int outVariable,
                          double coefficient) 
*/
int _init(PyObject* self, PyObject *args, PyObject *kwds) {
    PyCoefficient* py_coeff = reinterpret_cast<PyCoefficient*>(self);
    // failed to cast or self was not initialised - something horrible happened
    if (py_coeff == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                        "Failed to resolve self as PolynomialCoefficient in __init__");
        return -1;
    }
    // legal python to call initialised_object.__init__() to reinitialise, so
    // handle this case
    if (py_coeff->coeff != nullptr) {
        delete py_coeff->coeff;
        py_coeff->coeff = nullptr;
    }
    // read in arguments

    PyObject* py_index;
    int value_axis;
    double coefficient;
    static char *kwlist[] = {const_cast<char*>("index_by_vector"),
                             const_cast<char*>("output_axis"),
                             const_cast<char*>("coefficient_value"),
                             nullptr};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oid", kwlist,
                                    &py_index, &value_axis, &coefficient)) {
        return -1;
    }

    // convert from list to std::vector<int>
    // first check validity of coefficients
    if (!PyList_Check(py_index)) {
        PyErr_SetString(PyExc_TypeError,
                        "Failed to resolve index as a list");
        return -1;
    }
    size_t list_size = PyList_Size(py_index); // nb: size 0 is legal
    std::vector<int> index(list_size);
    // now loop over the rows
    for (size_t i = 0; i < list_size; ++i) {
        PyObject* py_value = PyList_GetItem(py_index, i);
        index[i] = int(PyLong_AsLong(py_value));
        if (PyErr_Occurred() != nullptr) { // not an int
            return -1;
        }
    }
    // now initialise the internal coeff
    try {
        py_coeff->coeff = new interpolation::PolynomialCoefficient(index, value_axis, coefficient);
    } catch (std::exception& exc) {
        PyErr_SetString(PyExc_RuntimeError, (&exc)->what());
        return -1;
    }
    return 0;
}

PyObject *_alloc(PyTypeObject *type, Py_ssize_t nitems) {
    void* void_coeff = malloc(sizeof(PyCoefficient));
    PyCoefficient* coeff = reinterpret_cast<PyCoefficient*>(void_coeff);
    coeff->coeff = nullptr;
    Py_REFCNT(coeff) = 1;
    Py_TYPE(coeff) = type;
    return reinterpret_cast<PyObject*>(coeff);
}

PyObject *_new(PyTypeObject *type, Py_ssize_t nitems) {
    return _alloc(type, nitems);
}

void _free(PyCoefficient * self) {
    if (self != nullptr) {
        if (self->coeff != nullptr)
            delete self->coeff;
        free(self);
    }
}

void _dealloc(PyCoefficient * self) {
    _free(self);
}

static PyMemberDef _members[] = {
{nullptr}
};

static PyMethodDef _methods[] = {
{nullptr}
};

std::string class_docstring =
std::string("PolynomialCoefficient docstring\n");

static PyTypeObject PyCoefficientType = {
    PyObject_HEAD_INIT(nullptr)
    "polynomial_coefficient.PolynomialCoefficient",         /*tp_name*/
    sizeof(PyCoefficient),           /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    class_docstring.c_str(),           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    _methods,           /* tp_methods */
    _members,           /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)_init,      /* tp_init */
    (allocfunc)_alloc,    /* tp_alloc, called by new */
    0, // (newfunc)_new,   /* tp_new */
    (freefunc)_free, /* tp_free, called by dealloc */
};

}  // namespace PyPolynomialCoefficient

const char* module_docstring =
        "polynomial_coefficient module contains the PolynomialCoefficient class";

static struct PyModuleDef polynomial_coefficient_def = {
    PyModuleDef_HEAD_INIT,
    "polynomial_coefficient",     /* m_name */
    module_docstring,  /* m_doc */
    -1,                  /* m_size */
    nullptr,    /* m_methods */
    nullptr,                /* m_reload */
    nullptr,                /* m_traverse */
    nullptr,                /* m_clear */
    nullptr,                /* m_free */
};

PyMODINIT_FUNC PyInit_polynomial_coefficient(void) {
    PyOpal::Globals::Initialise();
    PyPolynomialCoefficient::PyCoefficientType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PyPolynomialCoefficient::PyCoefficientType) < 0)
        return nullptr;

    PyObject* module = PyModule_Create(&polynomial_coefficient_def);
    if (module == nullptr)
        return nullptr;

    PyTypeObject* polynomial_coeff_type =
                          &PyPolynomialCoefficient::PyCoefficientType;
    Py_INCREF(polynomial_coeff_type);
    PyModule_AddObject(module, "PolynomialCoefficient",
                       reinterpret_cast<PyObject*>(polynomial_coeff_type));
    return module;
}


