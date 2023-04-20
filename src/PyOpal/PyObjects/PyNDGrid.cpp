//
// Python API for NDGrid
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
#include <structmember.h>

#include <exception>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "Utilities/OpalException.h"

#include "PyOpal/ExceptionTranslation.h"
#include "Classic/Fields/Interpolation/NDGrid.h"


namespace PyOpal {

namespace PyNDGrid {

using namespace interpolation;
namespace py = boost::python;

class Inform;
extern Inform *gmsg;

NDGrid initialiseVariableSpacing(boost::python::list gridCoordinates) {
    int gLength = boost::python::len(gridCoordinates);
    std::vector<std::vector<double> > coords(gLength);
    for (int i = 0; i < gLength; ++i) {
        int lineLength = boost::python::len(gridCoordinates[i]);
        coords[i] = std::vector<double>(lineLength);
        for (int j = 0; j < lineLength; ++j) {
            coords[i][j] = boost::python::extract<double>(gridCoordinates[i][j]);
        }
    }
    return NDGrid(coords);
}


NDGrid initialiseFixedSpacing(boost::python::list size,
                              boost::python::list spacing,
                              boost::python::list min) {
    int dim = boost::python::len(size);
    if (dim != boost::python::len(spacing)) {
        // error
    } else if (dim != boost::python::len(min)) {
        // error
    }
    std::vector<int> sizeVec(dim);
    std::vector<double> spacingVec(dim);
    std::vector<double> minVec(dim);
    for (int i = 0; i < dim; ++i) {
        sizeVec[i] = boost::python::extract<int>(size[i]);
        spacingVec[i] = boost::python::extract<double>(spacing[i]);
        minVec[i] = boost::python::extract<double>(min[i]);
    }

    return NDGrid(dim, &sizeVec[0], &spacingVec[0], &minVec[0]);
}

py::list coordVector(NDGrid& grid, int dimension) {
    if (dimension >= grid.getPositionDimension()) {
        throw OpalException("PyNDGrid::coordVector", 
                               "Dimension out of bounds");
    }
    std::vector<double> vec = grid.coordVector(dimension);
    py::list coord = py::list();
    for (size_t i = 0; i < vec.size(); ++i) {
        coord.append(vec[i]);
    }
    return coord;
}




const char* module_docstring = "ndgrid module for generating grids";

BOOST_PYTHON_MODULE(ndgrid) {
    ExceptionTranslation::registerExceptions();
    boost::python::class_<NDGrid>("NDGrid")
      .def("initialise_variable_spacing", &initialiseVariableSpacing)
      .staticmethod("initialise_variable_spacing")
      .def("initialise_fixed_spacing", &initialiseFixedSpacing)
      .staticmethod("initialise_fixed_spacing")
      .def("size", &NDGrid::size)
      .def("get_position_dimension", &NDGrid::getPositionDimension)
      .def("coord_vector", &coordVector)
    ;
}
} // namespace PyNDGrid
} // namespace PyOpal
