//
// Class NoMasterGraph
//   A simple empty master graph (no neighbors, every island works isolated).
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
#ifndef __NO_MASTER_GRAPH__
#define __NO_MASTER_GRAPH__

#include <set>

template < class TopoDiscoveryStrategy_t >
class NoMasterGraph : public TopoDiscoveryStrategy_t {

public:

    std::set<size_t> execute(size_t numMasters, size_t dimensions, size_t id,
                             int island_id) {
        std::set<size_t> empty;
        return empty;
    }
};

#endif
