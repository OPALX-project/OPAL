//
// Struct IndependentBitMutation
//   Mutate each gene with probability p
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

template <class T> struct IndependentBitMutation
{
    void mutate(boost::shared_ptr<T> ind, CmdArguments_t args) {

        const double probability =
            args->getArg<double>("gene-mutation-probability", 0.5);

        for(size_t i = 0; i < ind->genes_m.size(); i++) {
            double rval = static_cast<double>(rand() / (RAND_MAX + 1.0));
            if(rval < probability) {
                ind->new_gene(i);
            }
        }
    }
};