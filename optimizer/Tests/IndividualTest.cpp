//
// Test IndividualTest
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
#include "Optimizer/EA/Individual.h"
#include "Optimizer/EA/Individual.h"
#include "gtest/gtest.h"

#include "boost/smart_ptr.hpp"

#include <string>

namespace {

    // The fixture for testing class Foo.
    class IndividualTest : public ::testing::Test {
    protected:

        IndividualTest() {
            // You can do set-up work for each test here.
            lower_bound_ = -10.0;
            upper_bound_ = 19.79;
        }

        virtual ~IndividualTest() {
            // You can do clean-up work that doesn't throw exceptions here.
        }

        // If the constructor and destructor are not enough for setting up
        // and cleaning up each test, you can define the following methods:

        virtual void SetUp() {
            // Code here will be called immediately after the constructor (right
            // before each test).
        }

        virtual void TearDown() {
            // Code here will be called immediately after each test (right
            // before the destructor).
        }

        boost::shared_ptr<Individual> createIndividual(size_t num_genes, std::string constraint="") {

            Individual::bounds_t bounds;
            Individual::names_t names;
            Individual::constraints_t constraints;
            for(size_t i=0; i < num_genes; i++) {
                bounds.push_back(
                    std::pair<double, double>(lower_bound_, upper_bound_));
                names.push_back("dvar"+std::to_string(i));
            }
            if (constraint.empty() == false)
                constraints.insert(std::pair<std::string,Expressions::Expr_t*>
                                   ("constraint0",new Expressions::Expr_t(constraint)));

            boost::shared_ptr<Individual> ind(new Individual(bounds,names,constraints));
            return ind;
        }

        double lower_bound_;
        double upper_bound_;
    };

    TEST_F(IndividualTest, IndividualRespectsBounds) {

        size_t num_genes = 1;
        boost::shared_ptr<Individual> ind = createIndividual(num_genes);
        double gene = ind->genes_m[0];

        EXPECT_LE(lower_bound_, gene) << "gene should respect lower bound";
        EXPECT_GE(upper_bound_, gene) << "gene should respect upper bound";

        size_t my_size = ind->genes_m.size();
        EXPECT_EQ(static_cast<size_t>(num_genes), my_size)
            << "individual should only have " << num_genes << " gene(s)";

    }

    TEST_F(IndividualTest, IndividualHasCorrectNumberOfGenes) {

        size_t num_genes = 12;
        boost::shared_ptr<Individual> ind = createIndividual(num_genes);

        size_t my_size = ind->genes_m.size();
        EXPECT_EQ(static_cast<size_t>(num_genes), my_size)
            << "individual should only have " << num_genes << " gene(s)";

    }

    TEST_F(IndividualTest, IndividualRandomGene) {

        size_t num_genes = 1;
        boost::shared_ptr<Individual> ind = createIndividual(num_genes);
        double gene = ind->genes_m[0];
        double new_gene = ind->new_gene(0);

        EXPECT_NE(gene, new_gene) << "new gene should be different";
    }

    TEST_F(IndividualTest, IndividualConstraint) {

        size_t num_genes = 2;
        double half = (lower_bound_ + upper_bound_) / 2.;
        std::string constraint  = "(dvar0 + dvar1)/2. <=" + std::to_string(half);

        // create several individuals to test
        for (int i=0; i<10; i++) {
            boost::shared_ptr<Individual> ind = createIndividual(num_genes,constraint);
            double gene0 = ind->genes_m[0];
            double gene1 = ind->genes_m[1];
            EXPECT_LE((gene0+gene1)/2, half) << "constraint should be respected";
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}