#include <Python.h>
#include <structmember.h>

#include <boost/algorithm/string.hpp>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "Main.cpp"
//#include "PyOpal/Globals.h"
#include "mpi.h"
#include "Parser/Parser.h" // Classic
#include "OpalParser/OpalParser.h"
#include "AbstractObjects/OpalData.h"

std::string initialise_from_opal_file_docstring = 
std::string("Initialise from opal file\n")+
std::string("If you are getting an error message from openMPI, try\n")+
std::string("rebuilding the MPI library with --disable-dlopen switch\n");

extern "C" {
PyObject* initialise_from_opal_file(PyObject */*self*/, PyObject *args, PyObject *kwds) {
    static char *kwlist[] = {const_cast<char*>("file_name"),
                             NULL};
    char* value;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|s", kwlist,
                                    &value)) {
        return NULL;
    }
    const char* exe = "parser";
    char* argvr[3];
    // argv must be NULL terminated array (a week of my life figuring that one)
    argvr[0] = static_cast<char*>(malloc(sizeof('c')*7));
    argvr[1] = value;
    argvr[2] = NULL;
    strcpy(argvr[0], exe);
    //strcpy(argvr[1], value);
    try {
        opalMain(2, argvr);
    } catch (...) {
        std::string err = "Failed to initialise OPAL from file";
        PyErr_SetString(PyExc_ValueError, err.c_str());
    }
    Py_RETURN_NONE;
}
}

std::string list_objects_docstring = "List objects";

PyObject* list_objects(PyObject */*self*/, PyObject */*args*/, PyObject */*kwds*/) {
    OpalData::getInstance()->printNames(std::cout, "");

    Py_RETURN_NONE;
}

const char* module_docstring = "parser module parses the input";
 
static PyMethodDef _module_methods[] = {
{"initialise_from_opal_file", (PyCFunction)initialise_from_opal_file,
  METH_VARARGS|METH_KEYWORDS, initialise_from_opal_file_docstring.c_str()},
{"list_objects", (PyCFunction)list_objects,
  METH_VARARGS|METH_KEYWORDS, list_objects_docstring.c_str()},
{NULL, NULL, 0, NULL}
};

static struct PyModuleDef parserdef = {
    PyModuleDef_HEAD_INIT,
    "parser",     /* m_name */
    module_docstring,  /* m_doc */
    -1,                  /* m_size */
    _module_methods,    /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
};

PyMODINIT_FUNC PyInit_parser(void) {
    //PyOpal::Globals::Initialise();
    PyObject* module = PyModule_Create(&parserdef);
    return module;
}