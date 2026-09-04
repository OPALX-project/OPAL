//
// Template SolutionState_t
//   Text (de)serialization helpers for SolutionState_t containers exchanged
//   between MasterNode instances. Replaces boost::archive; requires that the
//   container's value_type provides writeState(std::ostream&) const and
//   readState(std::istream&).
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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
#ifndef __SOLUTION_STATE_IO__
#define __SOLUTION_STATE_IO__

#include <ostream>
#include <istream>

template <class SolutionState_t>
void writeSolutionState(std::ostream& os, const SolutionState_t& state) {
    os << state.size();
    for (const auto& ind : state) {
        os << ' ';
        ind.writeState(os);
    }
}

template <class SolutionState_t>
void readSolutionState(std::istream& is, SolutionState_t& state) {
    size_t n = 0;
    is >> n;
    state.resize(n);
    for (auto& ind : state) {
        ind.readState(is);
    }
}

#endif