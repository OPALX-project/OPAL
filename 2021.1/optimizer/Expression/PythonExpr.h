//
// Struct PythonExpression
//   Execute a python script using all given arguments (except for the first
//   being the script name) and return the value of the "result" variable in
//   the python script.
//   Expression arguments can be accessed in the Python script using the vector
//   variable "arguments", like i.e. in file "test.py":
//   \verbatim
//      result = 2.0 * arguments[0]
//    \endverbatim
//   to double the value passed as second argument to the Python expression.
//   For example the expression
//   \verbatim
//      EXPR="python("test.py", 15.0)";
//    \endverbatim
//   evaluates to a value of 30.0 (test.py contains the line shown above).
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef __PYTHON_EXPR_H__
#define __PYTHON_EXPR_H__

#include "boost/tuple/tuple.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "Util/Types.h"
#include "Util/OptPilotException.h"
#include "Expression/Expression.h"
#include "Expression/Parser/function.hpp"

using namespace boost::python;

struct PythonExpression {

    Expressions::Result_t operator()(client::function::arguments_t args) {
        std::string script = boost::get<std::string>(args[0]);
        std::vector<double> pargs;
        for(size_t i = 1; i < args.size(); i++) {
            pargs.push_back(boost::get<double>(args[i]));
        }

        try {
            Py_Initialize();
            object main_module    = import("__main__");
            object main_namespace = main_module.attr("__dict__");

            boost::python::class_<std::vector<double> >("PyVec")
                .def(boost::python::vector_indexing_suite<std::vector<double> >());
            main_namespace["arguments"] = pargs;

            object ignored = exec_file(script.c_str(), main_namespace);
            double res = extract<double>(main_namespace["result"]);

            Py_Finalize();
            return boost::make_tuple(res, true);

        } catch (error_already_set) {
            PyErr_Print();
            return boost::make_tuple(0.0, false);
        }
    }
};


#endif
