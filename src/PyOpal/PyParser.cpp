#include <Python.h>
#include <structmember.h>

#include <boost/python.hpp>
#include <boost/python/docstring_options.hpp>
#include <string>
#include <vector>

#include "Main.cpp"
#include "mpi.h"
#include "Parser/Parser.h" // Classic
#include "PyOpal/ExceptionTranslation.h"
//#include "PyOpal/Globals.h" // Not allowed: multiple definitions linker error
#include "PyOpal/PyParser.h"

namespace PyOpal {

namespace py = boost::python;

namespace parser {

py::object initialise_from_opal_file(std::string file_name) {
    char exe[] = {'p', 'a', 'r', 's', 'e', 'r', '\0'}; // surely not!
    std::vector<char> fname(file_name.size()+1);
    memcpy(&fname.front(), file_name.c_str(), file_name.size()+1); // urg
    char* argvr[3];
    // argv must be NULL terminated (a week of my life figuring that one)
    argvr[0] = exe;
    argvr[1] = fname.data();
    argvr[2] = NULL;
    strcpy(argvr[0], exe);
    opalMain(2, argvr);
    return py::object(); //PyNone
}

std::string list_objects(std::string regular_expression) {
    std::stringstream objectsString;
    OpalData::getInstance()->printNames(objectsString, regular_expression);
    return objectsString.str();
}

py::docstring_options options(true, false);

BOOST_PYTHON_MODULE(parser) {
    ExceptionTranslation::registerExceptions();
    py::scope().attr("__doc__") = parser_docstring.c_str();
    py::def("initialise_from_opal_file",
            initialise_from_opal_file,
            py::args("file_name"),
            initialise_from_opal_file_docstring.c_str()
    );
    py::def("list_objects",
            list_objects,
            list_objects_docstring.c_str()
    );
}
}
}
