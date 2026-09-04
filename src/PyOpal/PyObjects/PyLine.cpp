//
// Python API for the Line
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

#include <Python.h>

#include "Lines/Line.h"
#include "PyOpal/PyCore/ExceptionTranslation.h"
#include "PyOpal/PyCore/Globals.h"
#include "PyOpal/PyObjects/PyLine.h"

// using namespace boost::python;
namespace PyOpal {
    template <>
    std::vector<PyOpalObjectNS::AttributeDef>
        PyOpalObjectNS::PyOpalObject<TBeamline<Element> >::attributes = {
            {"L", "length", "", PyOpalObjectNS::DOUBLE},
            {"ORIGIN", "origin", "", PyOpalObjectNS::STRING},
            {"ORIENTATION", "orientation", "", PyOpalObjectNS::STRING},
            {"X", "x", "", PyOpalObjectNS::DOUBLE},
            {"Y", "y", "", PyOpalObjectNS::DOUBLE},
            {"Z", "z", "", PyOpalObjectNS::DOUBLE},
            {"THETA", "theta", "", PyOpalObjectNS::DOUBLE},
            {"PHI", "phi", "", PyOpalObjectNS::DOUBLE},
            {"PSI", "psi", "", PyOpalObjectNS::DOUBLE}};

    namespace PyLineNS {

        BOOST_PYTHON_MODULE(line) {
            ExceptionTranslation::registerExceptions();
            PyOpal::Globals::Initialise();
            PyLine aLine;
            auto lineClass = aLine.make_class("Line");
            // https://docs.python.org/3/library/collections.abc.html
            lineClass.def("__len__", &PyLine::getLength)
                .def("__getitem__", &PyLine::getElement)
                .def("__setitem__", &PyLine::setElement)
                //.def("__delitem__", &PyLine::removeElement)
                .def("append", &PyLine::append);
            lineClass.def("register", &PyLine::registerObject);
            lineClass.def("get_opal_name", &PyLine::getName);
            lineClass.def("set_opal_name", &PyLine::setName);
            aLine.addGetOpalElement(lineClass);

            // line is dependent on opal_element; all line elements are stored as
            // abstract opal_elements and we need boost to know how to do the
            // translation
            PyObject* mod = PyImport_ImportModule("pyopal.elements.opal_element");
            if (mod == nullptr) {
                PyErr_Print();
            }
        }

    }  // namespace PyLineNS
}  // namespace PyOpal
