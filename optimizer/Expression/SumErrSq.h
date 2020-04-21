//
// Struct SumErrSq
//  A simple expression computing the sum of all measurement errors (given as
//  first and third argument) for a variable (second argument) according to
//
//  \f[
//     result = \frac{1}{n} * \sqrt{\sum_{i=0}^n (measurement_i - value_i)^2}
//  \f]
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
#ifndef __SUMERRSQ_H__
#define __SUMERRSQ_H__

#include <map>
#include <string>
#include <fstream>
#include <iterator>

#include "boost/type_traits/remove_cv.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/SDDSReader.h"
#include "Expression/Parser/function.hpp"

class Measurement {
public:
    double spos;
    double measurement;

    friend std::istream & operator>>(std::istream & stream, Measurement & measurement);
};


struct SumErrSq {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {

        if (args.size() != 3) {
            throw OptPilotException("SumErrSq::operator()",
                                    "sumErrSq expects 3 arguments, " + std::to_string(args.size()) + " given");
        }

        std::string measurement_filename = boost::get<std::string>(args[0]);
        var_name_                        = boost::get<std::string>(args[1]);
        stat_filename_                   = boost::get<std::string>(args[2]);

        //FIXME: we could assume measurements don't change
        parseMeasurements(measurement_filename);
        bool is_valid = true;

        boost::scoped_ptr<SDDSReader> sim_stats(new SDDSReader(stat_filename_));
        try {
            sim_stats->parseFile();
        } catch (OptPilotException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
            is_valid = false;
        }

        double sum = 0;

        for(Measurement measurement : measurements_) {
            double sim_value = 0.0;
            try {
                sim_stats->getInterpolatedValue(
                        measurement.spos, var_name_, sim_value);
            } catch(OptPilotException &e) {
                std::cout << "Exception while getting value "
                          << "from SDDS file: " << e.what()
                          << std::endl;
                is_valid = false;
            }
            double val = measurement.measurement - sim_value;
            sum += val * val;
        }

        return boost::make_tuple(sqrt(sum/measurements_.size()), is_valid);
    }

private:

    std::vector<Measurement> measurements_;

    std::string var_name_;
    std::string stat_filename_;

    void parseMeasurements(std::string measurement_filename);

    // define a mapping to arguments in argument vector
    boost::tuple<std::string, std::string, std::string> argument_types;
};

#endif