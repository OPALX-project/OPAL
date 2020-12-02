#include <Python.h>
#include <structmember.h>

#include <string>
#include <exception>

#include <boost/python.hpp>

namespace PyOpal {
namespace py = boost::python;

namespace ExceptionTranslation {

void registerExceptions();

template <class T>
void translateException(T const& exception) {
    PyErr_SetString(PyExc_RuntimeError, exception.what());
}

template <class T>
void translateOpalException(T const& exception) {
    std::string msg = exception.what()+" in C++ method "+exception.where();
    PyErr_SetString(PyExc_RuntimeError, msg.c_str());
} 


}
}