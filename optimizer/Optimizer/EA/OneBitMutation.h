//
// Struct OneBitMutation
//   Mutate exactly one gene of an individual.
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

template <class T> struct OneBitMutation
{
    void mutate(boost::shared_ptr<T> ind, CmdArguments_t /*args*/) {

        int range = ind->genes_m.size();
        int position = static_cast<int>((rand() / (RAND_MAX + 1.0)) * range);
        ind->new_gene(position);
    }
};
