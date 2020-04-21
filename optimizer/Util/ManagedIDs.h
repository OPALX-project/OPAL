//
// Class ManagedIDs
//   Simple class to manage an ID pool.
//
//   Previously freed ID's are redistributed on following requests.
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
#include <cstddef>
#include <queue>

class ManagedIDs {

public:

    ManagedIDs() : next_free_(0)
    {}

    /// return next free ID
    size_t nextID() {

        size_t id = 0;

        if(freeids_.empty()) {
            id = next_free_;
            next_free_++;
        } else {
            id = freeids_.front();
            freeids_.pop();
        }

        return id;
    }


    /// free previously allocated ID
    void freeID(size_t id) {

        if(id == next_free_ - 1)
            next_free_--;
        else
            freeids_.push(id);
    }


private:

    /// queue to handle freed ID's
    std::queue<size_t> freeids_;

    /// next free ID
    size_t next_free_;

};
