#ifndef __NUMBER_OF_PEAKS_H__
#define __NUMBER_OF_PEAKS_H__

#include <string>

#include <cmath>

#include "boost/type_traits/remove_cv.hpp"
#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/PeakReader.h"
#include "Expression/Parser/function.hpp"

/**
 *  A simple expression to check the number of turns in a circular machine. It checks
 *  probe files (*.peaks) and counts the number of turns.
 */
struct NumberOfPeaks {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("NumberOfPeaks::operator()",
                                    "numberOfPeaks expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        sim_filename_  = boost::get<std::string>(args[0]);

        bool is_valid = true;

        boost::scoped_ptr<PeakReader> sim_peaks(new PeakReader(sim_filename_));
        std::size_t nPeaks = 0;

        try {
            sim_peaks->parseFile();
            nPeaks = sim_peaks->getNumberOfPeaks();
        } catch (OptPilotException &ex) {
            std::cout << "Caught exception: " << ex.what() << std::endl;
            is_valid = false;
        }


        return boost::make_tuple(nPeaks, is_valid);
    }

private:
    std::string sim_filename_;

    // define a mapping to arguments in argument vector
    boost::tuple<std::string> argument_types;
};

#endif
// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
