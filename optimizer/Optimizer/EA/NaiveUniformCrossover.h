#include "boost/smart_ptr.hpp"
#include "Util/CmdArguments.h"

/// decide for each gene if swapped with other gene
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
