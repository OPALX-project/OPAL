//
// Struct SeptumExpr
//   Objective to obtain a nice septum in cyclotron simulations.
//
// Copyright (c) 2019, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#ifndef __SEPTUM_EXPRESSION_H__
#define __SEPTUM_EXPRESSION_H__

#include <string>
#include <iostream>

#include "boost/variant/get.hpp"
#include "boost/variant/variant.hpp"
#include "boost/smart_ptr.hpp"

#include "Util/Types.h"
#include "Util/PeakReader.h"
#include "Expression/Parser/function.hpp"

#include "Util/ProbeHistReader.h"

struct SeptumExpr {

    static const std::string name;

    Expressions::Result_t operator()(client::function::arguments_t args) {
        if (args.size() != 1) {
            throw OptPilotException("SeptumExpr::operator()",
                                    "SeptumExpr expects 1 arguments, " + std::to_string(args.size()) + " given");
        }

        std::string probe = boost::get<std::string>(args[0]);

        bool is_valid = true;

        double result = 0.0;

        try {
            boost::scoped_ptr<PeakReader> sim_peaks(new PeakReader(probe + std::string(".peaks")));
            sim_peaks->parseFile();

            boost::scoped_ptr<ProbeHistReader> sim_hist(new ProbeHistReader(probe + std::string(".hist")));
            sim_hist->parseFile();

            double upperBound = 0.0;
            double lowerBound = 0.0;

            size_t nTurns = sim_peaks->getNumberOfPeaks();
            sim_peaks->getPeak(nTurns, upperBound);
            sim_peaks->getPeak(nTurns - 1, lowerBound);

            result = sim_hist->minimum(lowerBound, upperBound);

        } catch (OptPilotException &ex) {
            std::cout << "Exception while getting septum value "
                      << ex.what()
                      << std::endl;
            is_valid = false;
        }

        return boost::make_tuple(result, is_valid);
    }
};

#endif
