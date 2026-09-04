//
// Class Individual
//   Structure for an individual in the population holding genes and objective
//   values.
//
//   @see Types.h
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
#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "Expression/Expression.h"

class Individual {

public:
    /// representation of genes
    typedef std::vector<double> genes_t;
    /// gene names
    typedef std::vector<std::string> names_t;
    /// objectives array
    typedef std::vector<double> objectives_t;
    /// bounds on design variables
    typedef std::vector< std::pair<double, double> > bounds_t;
    /// constraints
    typedef Expressions::Named_t constraints_t;

    Individual()
    {}

    /// create a new individual and initialize with random genes
    Individual(bounds_t gene_bounds, names_t names, constraints_t constraints)
        : bounds_m(gene_bounds)
        , names_m(names)
        , constraints_m(constraints)
    {
        genes_m.resize(bounds_m.size(), 0.0);
        objectives_m.resize(bounds_m.size(), 0.0);

        // names should be equal length to bounds
        if (names_m.size() != bounds_m.size()) {
            // shouldn't happen
            std::cerr << "Individual::Individual(): names not equal length to bounds, shouldn't happen exiting" << std::endl;
            exit(1);
        }

        int iter = 0;
        while (iter < 100) { // somewhat arbitrary limit
            for (std::size_t i=0; i < bounds_m.size(); i++) {
                new_gene(i);
            }
            // check constraints
            bool allSatisfied = checkConstraints();
            if (allSatisfied == true) break;
            // else next try
            iter++;
        }
    }

    /// copy another individual
    Individual(std::shared_ptr<Individual> individual) {
        genes_m       =       genes_t(individual->genes_m);
        objectives_m  =  objectives_t(individual->objectives_m);
        bounds_m      =      bounds_t(individual->bounds_m);
        names_m       =       names_t(individual->names_m);
        constraints_m = constraints_t(individual->constraints_m);
        id_m          = individual->id_m;
    }

    /// write this individual's state to a stream, replaces boost::serialization
    void writeState(std::ostream& os) const {
        os << std::setprecision(std::numeric_limits<double>::max_digits10);
        os << id_m << ' ' << genes_m.size();
        for (double g : genes_m) os << ' ' << g;
        os << ' ' << objectives_m.size();
        for (double o : objectives_m) os << ' ' << o;
        os << '\n';
    }

    /// read this individual's state from a stream, replaces boost::serialization
    void readState(std::istream& is) {
        std::size_t ngenes = 0, nobjs = 0;
        is >> id_m >> ngenes;
        genes_m.resize(ngenes);
        for (double& g : genes_m) is >> g;
        is >> nobjs;
        objectives_m.resize(nobjs);
        for (double& o : objectives_m) is >> o;
    }

    /// initialize the gene with index gene_idx with a new random value
    /// contained in the specified gene boundaries
    double new_gene(std::size_t gene_idx) {
        double max = std::max(bounds_m[gene_idx].first, bounds_m[gene_idx].second);
        double min = std::min(bounds_m[gene_idx].first, bounds_m[gene_idx].second);
        double delta = std::abs(max - min);
        genes_m[gene_idx] = std::rand() / (RAND_MAX + 1.0) * delta + min;
        return genes_m[gene_idx];
    }
    /// test if individual within bounds and constraints
    bool viable() {
        return checkBounds() && checkConstraints();
    }

    /// genes of an individual
    genes_t      genes_m;
    /// values of objectives of an individual
    objectives_t objectives_m;
    /// id
    unsigned int id_m = 0;

private:
    /// check bounds
    bool checkBounds() {
        for (std::size_t i=0; i < bounds_m.size(); i++) {
             double value = genes_m[i];
             double max = std::max(bounds_m[i].first, bounds_m[i].second);
             double min = std::min(bounds_m[i].first, bounds_m[i].second);
             bool is_valid = (value >= min && value<= max);
             if (is_valid == false) {
                 return false;
             }
        }
        return true;
    }

    /// check if all constraints on design variables are checked
    bool checkConstraints() {
        for (auto namedConstraint : constraints_m) {
            Expressions::Expr_t *constraint = namedConstraint.second;
            std::set<std::string> req_vars = constraint->getReqVars();
            variableDictionary_t variable_dictionary;
            // fill variable dictionary. all required variables should be genes.
            for (std::string req_var : req_vars) {
                auto it = std::find(names_m.begin(),names_m.end(),req_var);
                if (it==names_m.end()) {
                    // should not happen
                    std::cerr << "Individual::checkConstraints(): " << req_var << " is not a design variable" << std::endl;
                    exit(1);
                }
                std::size_t gene_idx = std::distance(names_m.begin(),it);
                double value    = genes_m[gene_idx];
                variable_dictionary.insert(std::pair<std::string, double>(req_var, value));
            }
            Expressions::Result_t result =
                constraint->evaluate(variable_dictionary);

            double evaluation = std::get<0>(result);
            bool   is_valid   = std::get<1>(result);

            if (is_valid==false || evaluation==0) {
                return false;
            }
        }
        return true;
    }
    /// bounds on each gene
    bounds_t bounds_m;
    /// gene names
    names_t names_m;
    /// constraints that depend only on design variables
    constraints_t constraints_m;
};

#endif