//
// Struct NaiveUniformCrossover
//   Decide for each gene if swapped with other gene.
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

template <class T> struct NaiveUniformCrossover
{
    void crossover(boost::shared_ptr<T> ind1, boost::shared_ptr<T> ind2,
                   CmdArguments_t /*args*/) {

        Individual::genes_t genes_ind2 = ind2->genes_m;

        for(std::size_t i = 0; i < ind1->genes_m.size(); i++) {
            int choose = (int) (2.0 * (double) rand() / (RAND_MAX + 1.0));
            if(choose == 1) {
                ind2->genes_m[i] = ind1->genes_m[i];
                ind1->genes_m[i] = genes_ind2[i];
            }
        }
    }
};
