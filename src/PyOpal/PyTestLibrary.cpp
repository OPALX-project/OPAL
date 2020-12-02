#include <Python.h>
#include <structmember.h>

#include <string>
#include "PyOpal/PyTestLibrary.h"

std::string module_docstring =
                       "mylibrary module does nothing";

PyMODINIT_FUNC initmy_library(void) {
    PyObject* module = Py_InitModule3("my_library", NULL, module_docstring.c_str());
    if (module == NULL)
        return;
}

