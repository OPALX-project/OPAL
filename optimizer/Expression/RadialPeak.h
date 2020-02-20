#ifndef __RADIALPEAK_H__
#define __RADIALPEAK_H__

#include <string>

#include "boost/type_traits/remove_cv.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/PeakReader.h"
#include "Expression/Parser/function.hpp"

/**
 * A simple expression to get the n-th peak of a radial probe
 *
 */
struct RadialPeak {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 2) {
            throw OptPilotException("RadialPeak::operator()",
                                    "radialPeak expects 2 arguments, " + std::to_string(args.size()) + " given");
        }

        peak_filename_ = boost::get<std::string>(args[0]);
        turn_number_   = boost::get<double>(args[1]);

        bool is_valid = true;

        boost::scoped_ptr<PeakReader> sim_peaks(new PeakReader(peak_filename_));
        try {
            sim_peaks->parseFile();
        } catch (OptPilotException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
            is_valid = false;
        }

        double sim_radius = 0.0;
        try {
            sim_peaks->getPeak(turn_number_, sim_radius);
        } catch(OptPilotException &e) {
            std::cout << "Exception while getting value "
                      << "from peak file: " << e.what()
                      << std::endl;
            is_valid = false;
        }

        return boost::make_tuple(sim_radius, is_valid);
    }

private:

    std::string peak_filename_;
    int turn_number_;

    // define a mapping to arguments in argument vector
    boost::tuple<std::string, int> argument_types;
    // :FIXME: remove unused enum
#if 0
    enum {
          peak_filename
        , turn_number
    } argument_type_id;
#endif
};

#endif
// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
