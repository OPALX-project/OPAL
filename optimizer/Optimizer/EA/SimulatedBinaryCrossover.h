//
// Struct SimulatedBinaryCrossover
//   Deb (1995) Simulated Binary Crossover (SBX)
//   Respects interval schemata.
//   Offspring are symmetric around parent solutions.
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
#include <cmath>

#include "boost/smart_ptr.hpp"
#include "Util/CmdArguments.h"


template <class T> struct SimulatedBinaryCrossover
{
    void crossover(boost::shared_ptr<T> ind1, boost::shared_ptr<T> ind2,
                   CmdArguments_t args) {

        double nu_c = args->getArg<double>("simbin-crossover-nu", 2.0, false);

        for(std::size_t i = 0; i < ind1->genes_m.size(); i++) {

            double ui = (double) rand() / (RAND_MAX + 1.0);
            double beta_qi = 0.0;
            if(ui <= 0.5) {
                beta_qi = pow(2 * ui, 1.0/(nu_c + 1.0));
            } else {
                beta_qi = pow(1.0/(2 * (1.0 - ui)), 1.0/(nu_c + 1.0));
            }

            double ming = std::min(ind1->genes_m[i], ind2->genes_m[i]);
            double maxg = std::max(ind1->genes_m[i], ind2->genes_m[i]);

            ind1->genes_m[i] = 0.5 * ((1 + beta_qi) * ming + (1 - beta_qi) * maxg);
            ind2->genes_m[i] = 0.5 * ((1 - beta_qi) * ming + (1 + beta_qi) * maxg);
        }
    }
};
