//
// Namespace GlobalFunctions
//   Defines multiple expressions to be evaluated in
//   objectives.
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
#ifndef __GLOBAL_FUNCTIONS_H__
#define __GLOBAL_FUNCTIONS_H__

#include <string>
#include <cmath>

#include "boost/tuple/tuple.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"

#include "Util/Types.h"
#include "Util/SDDSReader.h"
#include "Util/OptPilotException.h"
#include "Expression/Expression.h"
#include "Expression/Parser/function.hpp"

namespace GlobalFunctions {

//FIXME
/*
struct _stat_sum {

    double operator()(client::function::arguments_t args) const {

        // get all arguments
        double start     = boost::get<double>(args[0]);
        double end       = boost::get<double>(args[1]);
        double step      = boost::get<double>(args[2]);
        std::string name = boost::get<std::string>(args[3]);

        boost::scoped_ptr<SDDSReader> sim_stats(new SDDSReader(stat_filename_));
        try {
            sim_stats->parseFile();
        } catch (OptPilotException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
        }

        double sum    = 0.0;
        bool is_valid = true;
        for(size_t i = static_cast<size_t>(start);
                   i < static_cast<size_t>(end); i++) {

            double sim_value;
            try {
                sim_stats->getInterpolatedValue(i, name, sim_value);
            } catch(OptPilotException &e) {
                std::cout << "Exception while getting value "
                          << "from SDDS file: " << e.what()
                          << std::endl;
                is_valid = false;
            }

            sum += sim_value;
        }

        return boost::make_tuple(value, is_valid);
    }
};
*/

struct _sqrt {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_sqrt::operator()",
                                    "sqrt expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(sqrt(value), true);
    }
};

struct _sq {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_sq::operator()",
                                    "sq expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(value * value, true);
    }
};

struct _pow {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 2) {
            throw OptPilotException("_pow::operator()",
                                    "pow expects 2 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        double exp   = boost::get<double>(args[1]);
        return boost::make_tuple(pow(value, exp), true);
    }
};

struct _exp {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_exp::operator()",
                                    "exp expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(exp(value), true);
    }
};

struct _log {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_log::operator()",
                                    "log expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(log(value), true);
    }
};



struct _ceil {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_ceil::operator()",
                                    "ceil expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(ceil(value), true);
    }
};

struct _fabs {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_fabs::operator()",
                                    "fabs expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(fabs(value), true);
    }
};

struct _floor {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_floor::operator()",
                                    "floor expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(floor(value), true);
    }
};

struct _fmod {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 2) {
            throw OptPilotException("_fmod::operator()",
                                    "fmod expects 2 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        double val2  = boost::get<double>(args[1]);
        return boost::make_tuple(fmod(value, val2), true);
    }
};



struct _sin {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_sin::operator()",
                                    "sin expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(sin(value), true);
    }
};

struct _asin {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_asin::operator()",
                                    "asin expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(asin(value), true);
    }
};

struct _cos {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_cos::operator()",
                                    "cos expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(cos(value), true);
    }
};

struct _acos {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_acos::operator()",
                                    "acos expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(acos(value), true);
    }
};

struct _tan {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_tan::operator()",
                                    "tan expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(tan(value), true);
    }
};

struct _atan {
    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("_atan::operator()",
                                    "atan expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        double value = boost::get<double>(args[0]);
        return boost::make_tuple(atan(value), true);
    }
};


static functionDictionary_t get() {

    typedef std::pair<std::string, client::function::type> funcEntry_t;
    functionDictionary_t funcs;

    client::function::type sqrt_; sqrt_ = _sqrt();
    funcs.insert(funcEntry_t("sqrt", sqrt_));
    client::function::type sq_; sq_ = _sq();
    funcs.insert(funcEntry_t("sq", sq_));
    client::function::type pow_; pow_ = _pow();
    funcs.insert(funcEntry_t("pow", pow_));
    client::function::type exp_; exp_ = _exp();
    funcs.insert(funcEntry_t("exp", exp_));
    client::function::type log_; log_ = _log();
    funcs.insert(funcEntry_t("log", log_));


    client::function::type ceil_; ceil_ = _ceil();
    funcs.insert(funcEntry_t("ceil", ceil_));
    client::function::type fabs_; fabs_ = _fabs();
    funcs.insert(funcEntry_t("fabs", fabs_));
    client::function::type floor_; floor_ = _floor();
    funcs.insert(funcEntry_t("floor", floor_));
    client::function::type fmod_; fmod_ = _fmod();
    funcs.insert(funcEntry_t("fmod", fmod_));


    client::function::type sin_; sin_ = _sin();
    funcs.insert(funcEntry_t("sin", sin_));
    client::function::type asin_; asin_ = _asin();
    funcs.insert(funcEntry_t("asin", asin_));
    client::function::type cos_; cos_ = _cos();
    funcs.insert(funcEntry_t("cos", cos_));
    client::function::type acos_; acos_ = _acos();
    funcs.insert(funcEntry_t("acos", acos_));
    client::function::type tan_; tan_ = _tan();
    funcs.insert(funcEntry_t("tan", tan_));
    client::function::type atan_; atan_ = _atan();
    funcs.insert(funcEntry_t("atan", atan_));

    return funcs;
}

}

#endif