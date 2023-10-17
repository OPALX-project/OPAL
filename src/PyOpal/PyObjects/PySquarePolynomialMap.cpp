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

#include <vector>
#include <string>

#include "Fields/Interpolation/MMatrix.h"
#include "Fields/Interpolation/SquarePolynomialVector.h"
#include "Fields/Interpolation/SolveFactory.h"
#include "Fields/Interpolation/LeastSquaresSolveFactory.h"

#include "PyOpal/Globals.h"
#include "PyOpal/PySquarePolynomialMap.h"
#include "PyOpal/PyPolynomialCoefficient.h"

namespace PySquarePolynomialMap {

using PyPolynomialCoefficient::PyCoefficient;
    
std::string get_coefficients_as_matrix_docstring =
std::string("Return the coefficients of the matrix. Takes no arguments\n\n")+
std::string("Return value is a list of lists, with polynomial coefficients\n")+
std::string("for y_i = sum_j a_{ij}*prod_n(x_n^j_n) forming the i, j term\n")+
std::string("in the matrix.\n");

PyObject* get_coefficients_as_matrix(PyObject *self, PyObject *args, PyObject *kwds) {
    PyPolynomialMap* py_map = reinterpret_cast<PyPolynomialMap*>(self);
    // failed to cast or self was not initialised - something horrible happened
    if (py_map == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                "Failed to resolve self as PolynomialMap");
        return nullptr;
    }
    if (py_map->map == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                "PolynomialMap not properly initialised");
        return nullptr;
    }
    interpolation::MMatrix<double> coefficients = py_map->map->GetCoefficientsAsMatrix();
    PyObject* py_coefficients = PyList_New(coefficients.num_row());
    // not safe from out of memory errors (etc)
    for (size_t i = 0; i < coefficients.num_row(); ++i) {
        PyObject* py_row = PyList_New(coefficients.num_col());
        for (size_t j = 0; j < coefficients.num_col(); ++j) {
            PyObject* py_value = PyFloat_FromDouble(coefficients(i+1, j+1));
            Py_INCREF(py_value);
            PyList_SetItem(py_row, j, py_value);
        }
        PyList_SetItem(py_coefficients, i, py_row);
        Py_INCREF(py_row);
    }
    Py_INCREF(py_coefficients);
    return py_coefficients;
}

std::string index_by_power_docstring =
std::string("Maps from the matrix index to the polynomial term.\n\n")+
std::string("  - col: index of the matrix column.\n")+
std::string("  - dim: dimension of the problem.\n")+
std::string("Return value is a list of ints of length dim. Terms in the\n")+
std::string("transfer matrix in the column given by col correspond to the\n")+
std::string("coefficient of the x_j1^k1 ... x_jn^kn term where kn are the\n")+
std::string("integers returned by this method.\n")+
std::string("For example, PolynomialMap.index_by_power(5, 2) returns\n")+
std::string("[2, 1] meaning the 7th column in a two-dimensional problem\n")+
std::string("represents the coefficients of x_0^2 x_1^1\n");


PyObject* index_by_power(PyObject *py_class, PyObject *args, PyObject *kwds) {
    PyTypeObject* py_map_type = reinterpret_cast<PyTypeObject*>(py_class);
    // failed to cast or self was not initialised - something horrible happened
    if (py_map_type == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                "Failed to resolve self as PolynomialMapType");
        return nullptr;
    }

    int col = 0;
    int dim = 0;
    static char *kwlist[] = {
        const_cast<char*>("col"),
        const_cast<char*>("dim"),
        nullptr
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ii", kwlist, &col, &dim)) {
        return nullptr;
    }
    if (col < 0) {
        PyErr_SetString(PyExc_ValueError, "col should be >= 0");
        return nullptr;
    }
    if (dim <= 0) {
        PyErr_SetString(PyExc_ValueError, "dim should be > 0");
        return nullptr;
    }
    std::vector<int> powers =
                  interpolation::SquarePolynomialVector::IndexByPower(col, dim);
    PyObject* py_powers = PyList_New(powers.size());
    Py_INCREF(py_powers);
    for (size_t i = 0; i < powers.size(); ++i) {
        PyObject* py_one_power = PyLong_FromSize_t(powers[i]);
        Py_INCREF(py_one_power);
        PyList_SetItem(py_powers, i, py_one_power);
    }
    return py_powers;
}

std::string evaluate_docstring =
std::string("Apply the mapping to a point.\n\n")+
std::string("  - point: list of floats with length equal to the point\n")+
std::string("    dimension of the mapping, corresponding to the abscissa.\n")+
std::string("Return value is a list of floats, with length equal to the\n")+
std::string("value dimension of the mapping; corresponding to the ordinates\n");

PyObject* evaluate(PyObject *self, PyObject *args, PyObject *kwds) {
    PyPolynomialMap* py_map = reinterpret_cast<PyPolynomialMap*>(self);
    // failed to cast or self was not initialised - something horrible happened
    if (py_map == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                "Failed to resolve self as PolynomialMap");
        return nullptr;
    }
    if (py_map->map == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                "PolynomialMap not properly initialised");
        return nullptr;
    }
    PyObject* py_point;
    static char *kwlist[] = {const_cast<char*>("point"), nullptr};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &py_point)) {
        return nullptr;
    }

    if (!PyList_Check(py_point)) {
        PyErr_SetString(PyExc_TypeError, "point was not a list");
        return nullptr;
    }
    if (PyList_Size(py_point) != py_map->map->PointDimension()) {
        PyErr_SetString(PyExc_TypeError, "point had wrong size");
        return nullptr;
    }
    std::vector<double> point(py_map->map->PointDimension());
    for (size_t i = 0; i < point.size(); ++i) {
        PyObject* point_i = PyList_GetItem(py_point, i);
        point[i] = PyFloat_AsDouble(point_i);
    }
    if (PyErr_Occurred()) // probably not a double in the list
        return nullptr;
    std::vector<double> value(py_map->map->ValueDimension());
    py_map->map->F(&point[0], &value[0]);
    PyObject* py_value = PyList_New(value.size());
    for (size_t i = 0; i < value.size(); ++i) {
        PyObject* value_i = PyFloat_FromDouble(value[i]);
        PyList_SetItem(py_value, i, value_i);
        Py_INCREF(value_i);
    }
    Py_INCREF(py_value);
    return py_value;
}

std::string exact_solve_docstring =
std::string("Find a mapping by solving for the matrix.\n\n")+
std::string("  - points: list, each entry containing a list of floats\n")+
std::string("  corresponding to abscissa, each with length PointDimension.\n")+
std::string("  The list should be as long as the number of coefficients in\n")+
std::string("  the polynomial.\n")+
std::string("  - values: list, each entry containing a list of floats\n")+
std::string("  corresponding to ordinates, each with length ValueDimension.\n")+
std::string("  The list should be as long as points.\n")+
std::string("  - polynomial_order: integer, >= 0, corresponding to the\n")+
std::string("  order of the fitted polynomial; 0 is constant, 1 is linear...\n")+
std::string("Returns a polynomial map.\n");

std::vector<std::vector<double> > get_vectors(PyObject* py_floats) {
    // convert from list of list of floats to std::vector< std::vector<double> >
    // first check validity of coefficients
     std::vector< std::vector<double> > data;
    if (!PyList_Check(py_floats)) {
        return std::vector<std::vector<double> >();
    }
    size_t num_rows = PyList_Size(py_floats);
    if (num_rows == 0) {
        return std::vector<std::vector<double> >();
    }
    size_t num_cols = 0;
    // now loop over the rows
    for (size_t i = 0; i < num_rows; ++i) {
        PyObject* row = PyList_GetItem(py_floats, i);
        if (!PyList_Check(row)) {
            return std::vector<std::vector<double> >();
        }
        // do some initialisation in the first row
        if (i == 0) {
            num_cols = PyList_Size(row);
            if (num_cols == 0) {
                return std::vector<std::vector<double> >();
            }
            data = std::vector< std::vector<double> >(num_rows, std::vector<double>(num_cols));
        }
        // now loop over columns
        if (PyList_Size(row) != static_cast<int>(num_cols)) {
            return std::vector<std::vector<double> >();
        }
        for (size_t j = 0; j < num_cols; ++j) {
            PyObject* py_value = PyList_GetItem(row, j);
            data.at(i).at(j) = PyFloat_AsDouble(py_value);
            if (PyErr_Occurred() != nullptr) // not a float
                return std::vector<std::vector<double> >();
        }
    }
    return data;
}

PyObject* exact_solve(PyObject *py_class, PyObject *args, PyObject *kwds) {
    PyTypeObject* py_map_type = reinterpret_cast<PyTypeObject*>(py_class);
    // failed to cast or self was not initialised - something horrible happened
    if (py_map_type == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                "Failed to resolve self as PolynomialMapType");
        return nullptr;
    }

    PyObject* py_points = nullptr;
    PyObject* py_values = nullptr;
    int polynomial_order = 0;
    PyObject* py_error_matrix = Py_None; // borrowed reference
    static char *kwlist[] = {
        const_cast<char*>("points"),
        const_cast<char*>("values"),
        const_cast<char*>("polynomial_order"),
        const_cast<char*>("error_matrix"),
        nullptr
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOi|O", kwlist, &py_points,
            &py_values, &polynomial_order, &py_error_matrix)) {
        return nullptr;
    }
    if (polynomial_order < 0) {
        PyErr_SetString(PyExc_ValueError, "polynomial_order should be >= 0");
        return nullptr;
    }
    std::vector<std::vector<double> > points = get_vectors(py_points);
    if (points.size() == 0) {
        PyErr_SetString(PyExc_ValueError, "Failed to evaluate points");
        return nullptr;
    }
    std::vector<std::vector<double> > values = get_vectors(py_values);
    if (values.size() == 0) {
        PyErr_SetString(PyExc_ValueError, "Failed to evaluate values");
        return nullptr;
    }
    if (points.size() != values.size()) {
        PyErr_SetString(PyExc_ValueError, "points misaligned with values");
        return nullptr;
    }
    interpolation::SquarePolynomialVector* test_map = nullptr;
    try {
        std::vector<std::vector<double> > no_derivs;
        std::vector<std::vector<int> > no_indices;
        size_t dim = points[0].size();
        interpolation::SolveFactory solve(polynomial_order, polynomial_order, 
                                          dim, dim, points, no_derivs, no_indices);
        test_map = solve.PolynomialSolve(values, no_derivs);
    } catch (GeneralClassicException& exc) {
        PyErr_SetString(PyExc_ValueError, exc.what().c_str());
        return nullptr;
    }
    PyObject* py_map_obj = _alloc(py_map_type, 0);
    PyPolynomialMap* py_map = reinterpret_cast<PyPolynomialMap*>(py_map_obj);
    py_map->map = test_map;
    return py_map_obj;
}

std::string least_squares_docstring =
std::string("Find a mapping by solving for the matrix.\n\n")+
std::string("  - points: list, each entry containing a list of floats\n")+
std::string("  corresponding to abscissa, each with length PointDimension.\n")+
std::string("  The list should be at least as long as the number of \n")+
std::string("  coefficients in the polynomial.\n")+
std::string("  - values: list, each entry containing a list of floats\n")+
std::string("  corresponding to ordinates, each with length ValueDimension.\n")+
std::string("  The list should be as long as points.\n")+
std::string("  - polynomial_order: integer, >= 0, corresponding to the\n")+
std::string("  order of the fitted polynomial; 0 is constant, 1 is linear...\n")+
std::string("  - coefficients: list of PolynomialCoefficients. Fix these\n")+
std::string("    coefficients to some \n")+
std::string("  - weights: list, of length = 0, corresponding to\n")+
std::string("Returns a polynomial map.\n");


PyObject* least_squares(PyObject *py_class, PyObject *args, PyObject *kwds) {
    PyTypeObject* py_map_type = reinterpret_cast<PyTypeObject*>(py_class);
    // failed to cast or self was not initialised - something horrible happened
    if (py_map_type == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                "Failed to resolve self as PolynomialMapType");
        return nullptr;
    }

    PyObject* py_points = nullptr;
    PyObject* py_values = nullptr;
    int polynomial_order = 0;
    PyObject* py_coefficients = nullptr;
    PyObject* py_weights = nullptr;
    static char *kwlist[] = {
        const_cast<char*>("points"),
        const_cast<char*>("values"),
        const_cast<char*>("polynomial_order"),
        const_cast<char*>("polynomial_coefficients"),
        const_cast<char*>("weights"),
        nullptr
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOi|OO", kwlist, &py_points,
            &py_values, &polynomial_order, &py_coefficients, &py_weights)) {
        return nullptr;
    }
    if (polynomial_order < 0) {
        PyErr_SetString(PyExc_ValueError, "polynomial_order should be >= 0");
        return nullptr;
    }
    std::vector<std::vector<double> > points = get_vectors(py_points);
    std::vector<std::vector<double> > values = get_vectors(py_values);
    std::vector<interpolation::PolynomialCoefficient> coeff;
    std::vector<double> weights;
    // weights
    if (py_weights != nullptr) {
        if (!PyList_Check(py_weights)) {
            PyErr_SetString(PyExc_TypeError,
                            "Failed to resolve weights as a list");
            return nullptr;
        }
        size_t list_size = PyList_Size(py_weights); // nb: size 0 is legal
        weights = std::vector<double>(list_size);
        for (size_t i = 0; i < list_size; ++i) {
            PyObject* py_value = PyList_GetItem(py_weights, i);
            weights[i] = int(PyFloat_AsDouble(py_value));
            if (PyErr_Occurred() != nullptr) { // not an int
                return nullptr;
            }
        }
    }
    // coefficients
    if (py_coefficients != nullptr) {
        if (!PyList_Check(py_coefficients)) {
            PyErr_SetString(PyExc_TypeError,
                            "Failed to resolve coefficients as a list");
            return nullptr;
        }
        size_t list_size = PyList_Size(py_coefficients); // nb: size 0 is legal
        for (size_t i = 0; i < list_size; ++i) {
            PyObject* py_value = PyList_GetItem(py_coefficients, i);
            PyCoefficient* py_coeff = reinterpret_cast<PyCoefficient*>(py_value);
            if (py_coeff == nullptr) {
                PyErr_SetString(PyExc_TypeError,
                    "Failed to resolve list item as a PolynomialCoefficient.");
                return nullptr;
            }
            coeff.push_back(*(py_coeff->coeff));
        }

    }
    interpolation::SquarePolynomialVector* test_map = nullptr;
    try {
        interpolation::LeastSquaresSolveFactory solver(polynomial_order, points);
        solver.setWeights(weights);
        solver.setCoefficients(coeff);
        test_map = new interpolation::SquarePolynomialVector(solver.solve(values));
    } catch (GeneralClassicException& exc) {
        PyErr_SetString(PyExc_ValueError, exc.what().c_str());
        return nullptr;
    }
    PyObject* py_map_obj = _alloc(py_map_type, 0);
    PyPolynomialMap* py_map = reinterpret_cast<PyPolynomialMap*>(py_map_obj);
    py_map->map = test_map;
    return py_map_obj;
}



int _init(PyObject* self, PyObject *args, PyObject *kwds) {
    PyPolynomialMap* py_map = reinterpret_cast<PyPolynomialMap*>(self);
    // failed to cast or self was not initialised - something horrible happened
    if (py_map == nullptr) {
        PyErr_SetString(PyExc_TypeError,
                        "Failed to resolve self as PolynomialMap in __init__");
        return -1;
    }
    // legal python to call initialised_object.__init__() to reinitialise, so
    // handle this case
    if (py_map->map != nullptr) {
        delete py_map->map;
        py_map->map = nullptr;
    }
    // read in arguments
    int point_dim;
    PyObject* py_coefficients;
    static char *kwlist[] = {const_cast<char*>("point_dimension"),
                             const_cast<char*>("coefficients"), nullptr};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iO", kwlist,
                                         &point_dim, &py_coefficients)) {
        return -1;
    }

    // convert from list of list of floats to Matrix<double>
    // first check validity of coefficients
    if (!PyList_Check(py_coefficients)) {
        PyErr_SetString(PyExc_TypeError,
                        "Failed to resolve coefficients as a list");
        return -1;
    }
    size_t num_rows = PyList_Size(py_coefficients);
    if (num_rows == 0) {
        PyErr_SetString(PyExc_ValueError, "coefficients was empty");
        return -1;
    }
    size_t num_cols = 0;
    interpolation::MMatrix<double> coefficients;
    // now loop over the rows
    for (size_t i = 0; i < num_rows; ++i) {
        PyObject* row = PyList_GetItem(py_coefficients, i);
        if (!PyList_Check(py_coefficients)) {
            PyErr_SetString(PyExc_TypeError,
                            "Failed to resolve coefficients row as a list");
            return -1;
        }
        // do some initialisation in the first row
        if (i == 0) {
            num_cols = PyList_Size(row);
            if (num_cols == 0) {
                PyErr_SetString(PyExc_ValueError, "coefficients row was empty");
                return -1;
            }
            coefficients = interpolation::MMatrix<double>(num_rows, num_cols);
        }
        // now loop over columns
        if (PyList_Size(row) != static_cast<int>(num_cols)) {
                PyErr_SetString(PyExc_ValueError,
                                "coefficients row had wrong number of elements");
        }
        for (size_t j = 0; j < num_cols; ++j) {
            PyObject* py_value = PyList_GetItem(row, j);
            coefficients(i+1, j+1) = PyFloat_AsDouble(py_value);
            if (PyErr_Occurred() != nullptr) // not a float
                return -1;
        }
    }

    // now initialise the internal map
    try {
        py_map->map = new interpolation::SquarePolynomialVector(point_dim, coefficients);
    } catch (std::exception& exc) {
        PyErr_SetString(PyExc_RuntimeError, (&exc)->what());
        return -1;
    }
    return 0;
}

PyObject *_alloc(PyTypeObject *type, Py_ssize_t nitems) {
    void* void_map = malloc(sizeof(PyPolynomialMap));
    PyPolynomialMap* map = reinterpret_cast<PyPolynomialMap*>(void_map);
    map->map = nullptr;
    Py_REFCNT(map) = 1;
    Py_TYPE(map) = type;
    return reinterpret_cast<PyObject*>(map);
}

PyObject *_new(PyTypeObject *type, Py_ssize_t nitems) {
    return _alloc(type, nitems);
}

void _dealloc(PyPolynomialMap * self) {
    _free(self);
}

void _free(PyPolynomialMap * self) {
    if (self != nullptr) {
        if (self->map != nullptr)
            delete self->map;
        free(self);
    }
}

static PyMemberDef _members[] = {
{nullptr}
};

static PyMethodDef _methods[] = {
{"get_coefficients_as_matrix", (PyCFunction)get_coefficients_as_matrix,
  METH_VARARGS|METH_KEYWORDS, get_coefficients_as_matrix_docstring.c_str()},
{"evaluate", (PyCFunction)evaluate,
  METH_VARARGS|METH_KEYWORDS, evaluate_docstring.c_str()},
{"exact_solve", (PyCFunction)exact_solve,
  METH_CLASS|METH_VARARGS|METH_KEYWORDS, exact_solve_docstring.c_str()},
{"least_squares", (PyCFunction)least_squares,
  METH_CLASS|METH_VARARGS|METH_KEYWORDS, least_squares_docstring.c_str()},
{"index_by_power", (PyCFunction)index_by_power,
  METH_CLASS|METH_VARARGS|METH_KEYWORDS, index_by_power_docstring.c_str()},
{nullptr}
};

std::string class_docstring =
std::string("PolynomialMap provides routines to calculate multivariate \n")+
std::string("polynomials.\n\n")+
std::string("__init__()\n")+
std::string("    Takes two arguments.\n")+
std::string("    - point_dim: integer which defines the dimension of the\n")+
std::string("      points (abscissa)\n")+
std::string("    - coefficients: list of lists of floats which define the\n")+
std::string("      polynomial\n")+
std::string("The value dimension of the PolynomialMap is the number of rows\n")+
std::string("coefficients matrix\n");

static PyTypeObject PyPolynomialMapType = {
    PyObject_HEAD_INIT(nullptr)
    "polynomial_map.PolynomialMap",         /*tp_name*/
    sizeof(PyPolynomialMap),           /*tp_basicsize*/
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
    _methods,                  /* tp_methods */
    _members,                  /* tp_members */
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

}  // namespace PyPolynomialMap

const char* module_docstring =
                       "polynomial_map module contains the PolynomialMap class";

static struct PyModuleDef polynomial_map_def = {
    PyModuleDef_HEAD_INIT,
    "polynomial_map",     /* m_name */
    module_docstring,  /* m_doc */
    -1,                  /* m_size */
    nullptr,    /* m_methods */
    nullptr,                /* m_reload */
    nullptr,                /* m_traverse */
    nullptr,                /* m_clear */
    nullptr,                /* m_free */
};

PyMODINIT_FUNC PyInit_polynomial_map(void) {
    PyOpal::Globals::Initialise();
    PySquarePolynomialMap::PyPolynomialMapType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PySquarePolynomialMap::PyPolynomialMapType) < 0)
        return nullptr;

    PyObject* module = PyModule_Create(&polynomial_map_def);
    if (module == nullptr)
        return nullptr;

    PyTypeObject* polynomial_map_type = &PySquarePolynomialMap::PyPolynomialMapType;
    Py_INCREF(polynomial_map_type);
    PyModule_AddObject(module, "PolynomialMap",
                       reinterpret_cast<PyObject*>(polynomial_map_type));

    return module;
}


