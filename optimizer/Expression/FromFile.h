#ifndef __FROMFILE_H__
#define __FROMFILE_H__

#include <string>
#include <map>
#include <set>
#include <fstream>
#include <iterator>

#include "boost/tuple/tuple.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/OptPilotException.h"
#include "Expression/Parser/function.hpp"


/**
 *  Simple functor that reads vector data from a file. If the file contains
 *  more than one value the sum is returned.
 *  \f[
 *    result = \sum_{i=0}^n value_i
 *  \f]
 */
struct FromFile {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("FromFile::operator()",
                                    "fromFile expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        filename_   = boost::get<std::string>(args[0]);

        double sum = 0;
        bool is_valid = true;

        try {
            readValues();
        } catch(const OptPilotException& e) {
            return boost::make_tuple(0.0, false);
        }

        for(double obj_value : values_)
            sum += obj_value;

        return boost::make_tuple(sum, is_valid);
    }


private:

    std::vector<double> values_;

    std::string filename_;

    void readValues();

};

#endif