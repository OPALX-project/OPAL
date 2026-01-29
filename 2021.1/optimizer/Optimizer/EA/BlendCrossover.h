//
// Struct BlendCrossover
//   BLX-alpha (interval schemata)
//   Eshelman and Schaffer (1993)
//   Pick random solution in interval
//
//     [ x_i^(1,t) - \alpha(x_i^(2,t) - x_i^(1,t)),
//       x_i^(2,t) + \alpha((x_i^(2,t) - x_i^(1,t)) ]
//
//   at generation t.
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
#include "boost/smart_ptr.hpp"
#include "Util/CmdArguments.h"
#include <cmath>


template <class T> struct BlendCrossover
{
    void crossover(boost::shared_ptr<T> ind1, boost::shared_ptr<T> ind2,
                   CmdArguments_t /*args*/) {

        // BLX-0.5 performs better than BLX operators with any other \alpha
        // value
        const double alpha = 0.5;

        for(size_t i = 0; i < ind1->genes_m.size(); i++) {

            double ming = std::min(ind1->genes_m[i], ind2->genes_m[i]);
            double maxg = std::max(ind1->genes_m[i], ind2->genes_m[i]);
            double gamma1 = (1 + 2 * alpha) *
                static_cast<double>(rand() / (RAND_MAX + 1.0)) - alpha;
            double gamma2 = (1 + 2 * alpha) *
                static_cast<double>(rand() / (RAND_MAX + 1.0)) - alpha;
            ind1->genes_m[i] = (1 - gamma1) * ming + gamma1 * maxg;
            ind2->genes_m[i] = (1 - gamma2) * ming + gamma2 * maxg;
        }
    }
};

