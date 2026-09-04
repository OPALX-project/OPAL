//
// Class SampleIndividual
//   Structure for an individual in the population holding genes values.
//
// Copyright (c) 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
//                     Yves Ineichen, ETH Zürich
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
#ifndef __SAMPLE_INDIVIDUAL_H__
#define __SAMPLE_INDIVIDUAL_H__

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include "Utilities/OpalException.h"

class SampleIndividual {

public:
    /// representation of genes
    typedef std::vector<double> genes_t;
    /// gene names
    typedef std::vector<std::string> names_t;
    /// objectives array
    typedef std::vector<double> objectives_t;

    SampleIndividual()
    {}

    SampleIndividual(names_t names)
        : names_m(names)
    {
        genes.resize(names.size(), 0.0);
    }

    /// write this individual's state to a stream, replaces boost::serialization
    void writeState(std::ostream& os) const {
        os << std::setprecision(std::numeric_limits<double>::max_digits10);
        os << id << ' ' << genes.size();
        for (double g : genes) os << ' ' << g;
        os << ' ' << objectives.size();
        for (double o : objectives) os << ' ' << o;
        os << '\n';
    }

    /// read this individual's state from a stream, replaces boost::serialization
    void readState(std::istream& is) {
        std::size_t ngenes = 0, nobjs = 0;
        is >> id >> ngenes;
        genes.resize(ngenes);
        for (double& g : genes) is >> g;
        is >> nobjs;
        objectives.resize(nobjs);
        for (double& o : objectives) is >> o;
    }

    /// genes of an individual
    genes_t      genes;
    /// values of objectives of an individual
    objectives_t objectives;
    /// id
    unsigned int id = 0;

    int getIndex(std::string name) {
        auto res = std::find(std::begin(names_m), std::end(names_m), name);

        if (res == std::end(names_m)) {
            throw OpalException("SampleIndividual::getIndex()",
                                "Variable '" + name + "' not contained.");
        }
        return std::distance(std::begin(names_m), res);
    }

    std::string getName(std::size_t i) {
        return names_m[i];
    }

    void print(std::ostream& out) const {
        out << std::setw(8) << id << std::endl;
        for (unsigned int i = 0; i < genes.size(); ++ i) {
            out << names_m[i] << ": " << genes[i] << std::endl;
        }
    }
private:
    /// gene names
    names_t names_m;
};

inline
std::ostream& operator<<(std::ostream& out, const SampleIndividual& ind) {
    ind.print(out);

    return out;
}
#endif