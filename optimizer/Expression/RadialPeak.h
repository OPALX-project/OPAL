#ifndef __RADIALPEAK_H__
#define __RADIALPEAK_H__

#include <string>
#include <tuple>
#include <variant>

#include "Util/Types.h"
#include "Util/PeakReader.h"
#include "Expression/Parser/function.hpp"

/**
 * A simple expression to get the n-th peak of a radial probe
 *
 */
struct RadialPeak {

    static const std::string name;

    RadialPeak()
        : argument_types("peak_filename", 0)
    {}

    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 2) {
            throw OptPilotException("RadialPeak::operator()",
                                    "radialPeak expects 2 arguments, " + std::to_string(args.size()) + " given");
        }

        peak_filename_ = std::get<std::string>(args[0]);
        turn_number_   = static_cast<int>(std::get<double>(args[1]));

        bool is_valid = true;
        double sim_radius = 0.0;

        try {
            PeakReader sim_peaks(peak_filename_);
            sim_peaks.parseFile();
            sim_peaks.getPeak(turn_number_, sim_radius);
        } catch (const OptPilotException& e) {
            std::cerr << "RadialPeak exception: " << e.what() << std::endl;
            is_valid = false;
        }

        return std::make_tuple(sim_radius, is_valid);
    }

private:

    std::string peak_filename_;
    int turn_number_;

    // define a mapping to arguments in argument vector
    std::tuple<std::string, int> argument_types;
};

#endif
