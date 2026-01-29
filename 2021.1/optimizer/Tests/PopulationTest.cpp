//
// Test PopulationTest
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
#include "Optimizer/EA/Population.h"
#include "Optimizer/EA/Individual.h"
#include "gtest/gtest.h"

#include "boost/smart_ptr.hpp"

namespace {

    // The fixture for testing class Foo.
    class PopulationTest : public ::testing::Test {
    protected:

        PopulationTest() {
            // You can do set-up work for each test here.
            population_.reset(new Population<Individual>());
        }

        virtual ~PopulationTest() {
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
            // population_->clean_population();
        }

        boost::shared_ptr<Individual> createIndividual(size_t num_genes) {

            Individual::bounds_t bounds;
            Individual::names_t names;
            Individual::constraints_t constraints;
            for(size_t i=0; i < num_genes; i++) {
                bounds.push_back(std::pair<double, double>(0.0, 10.0));
                names.push_back("dvar"+std::to_string(i));
            }

            boost::shared_ptr<Individual> ind(new Individual(bounds,names,constraints));
            return ind;
        }

        // Objects declared here can be used by all tests in the test case
        boost::scoped_ptr< Population<Individual> > population_;
    };

    TEST_F(PopulationTest, AddOneIndividual) {

        boost::shared_ptr<Individual> ind = createIndividual(1);
        unsigned int id = population_->add_individual(ind);
        double gene = ind->genes_m[0];
        double obj  = ind->objectives_m[0];

        EXPECT_EQ(static_cast<size_t>(0), id) << "first individuals id should be 0";

        boost::shared_ptr<Individual> tmp = population_->get_individual(id);
        EXPECT_EQ(0, tmp.get()) << "no committed individuals after insert";

        tmp = population_->get_staging(id);
        EXPECT_EQ(gene, tmp->genes_m[0])      << "gene should have specified value";
        EXPECT_EQ(obj,  tmp->objectives_m[0]) << "objective should have specified value";

        size_t my_size = population_->size();
        EXPECT_EQ(static_cast<size_t>(0), my_size)
            << "no committed individuals, population should still be 0";

    }

    TEST_F(PopulationTest, CommitOneIndividual) {

        boost::shared_ptr<Individual> ind = createIndividual(1);
        unsigned int id = population_->add_individual(ind);
        double gene = ind->genes_m[0];
        double obj  = ind->objectives_m[0];

        EXPECT_EQ(static_cast<size_t>(0), id) << "first individuals id should be 0";

        population_->commit_individuals();

        boost::shared_ptr<Individual> tmp = population_->get_staging(id);
        EXPECT_EQ(0, tmp.get()) << "no staging individuals after commit";

        tmp = population_->get_individual(id);
        EXPECT_EQ(gene, tmp->genes_m[0]);
        EXPECT_EQ(obj,  tmp->objectives_m[0]);

        size_t my_size = population_->size();
        EXPECT_EQ(static_cast<size_t>(1), my_size);

    }

    TEST_F(PopulationTest, KeepIndividuals) {

        boost::shared_ptr<Individual> ind1 = createIndividual(1);
        boost::shared_ptr<Individual> ind2 = createIndividual(1);
        boost::shared_ptr<Individual> ind3 = createIndividual(1);
        boost::shared_ptr<Individual> ind4 = createIndividual(1);

        size_t id0 = population_->add_individual(ind1);
        EXPECT_EQ(static_cast<size_t>(0), id0);
        size_t id1 = population_->add_individual(ind2);
        EXPECT_EQ(static_cast<size_t>(1), id1);
        size_t id2 = population_->add_individual(ind3);
        EXPECT_EQ(static_cast<size_t>(2), id2);
        size_t id3 = population_->add_individual(ind4);
        EXPECT_EQ(static_cast<size_t>(3), id3);

        population_->commit_individuals();

        std::set<unsigned int> survivors;
        survivors.insert(id1);
        survivors.insert(id3);

        population_->keepSurvivors(survivors);

        size_t size = population_->size();
        EXPECT_EQ(survivors.size(), size);

        boost::shared_ptr<Individual> tmp = population_->get_individual(id1);
        EXPECT_EQ(ind2->genes_m[0], tmp->genes_m[0]);
        EXPECT_EQ(ind2->objectives_m[0], tmp->objectives_m[0]);
    }

    TEST_F(PopulationTest, IDsContinuous) {

        boost::shared_ptr<Individual> ind0 = createIndividual(1);
        boost::shared_ptr<Individual> ind1 = createIndividual(1);
        boost::shared_ptr<Individual> ind2 = createIndividual(1);
        boost::shared_ptr<Individual> ind3 = createIndividual(1);

        size_t id0 = population_->add_individual(ind0);
        EXPECT_EQ(static_cast<size_t>(0), id0);
        size_t id1 = population_->add_individual(ind1);
        EXPECT_EQ(static_cast<size_t>(1), id1);
        size_t id2 = population_->add_individual(ind2);
        EXPECT_EQ(static_cast<size_t>(2), id2);
        size_t id3 = population_->add_individual(ind3);
        EXPECT_EQ(static_cast<size_t>(3), id3);

        unsigned int individual_to_be_removed_id = 1;
        population_->remove_individual(ind1);

        boost::shared_ptr<Individual> newind = createIndividual(1);
        unsigned int id_new = population_->add_individual(newind);
        EXPECT_EQ(individual_to_be_removed_id, id_new);

        boost::shared_ptr<Individual> tmp = population_->get_staging(id_new);
        EXPECT_EQ(newind->genes_m[0], tmp->genes_m[0]);
        EXPECT_EQ(newind->objectives_m[0], tmp->objectives_m[0]);
    }

    TEST_F(PopulationTest, FindNonExistingStaging) {

        boost::shared_ptr<Individual> tmp = population_->get_staging(124);
        EXPECT_EQ(0, tmp.get());
    }

    TEST_F(PopulationTest, FindNonExistingIndividual) {

        boost::shared_ptr<Individual> tmp = population_->get_individual(124);
        EXPECT_EQ(0, tmp.get());
    }

    TEST_F(PopulationTest, RepresentedCheck) {

        boost::shared_ptr<Individual> ind = createIndividual(1);
        population_->add_individual(ind);

        population_->add_individual(ind);
        population_->commit_individuals();

        bool represented = population_->isRepresentedInPopulation(ind->genes_m);
        EXPECT_TRUE(represented);

        std::vector<double> tmp_genes;
        tmp_genes.push_back(12312);
        represented = population_->isRepresentedInPopulation(tmp_genes);
        EXPECT_FALSE(represented);
    }

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}