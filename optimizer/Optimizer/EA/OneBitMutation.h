#include "boost/smart_ptr.hpp"
#include "Util/CmdArguments.h"

/// Mutate exactly one gene of an individual.
template <class T> struct OneBitMutation
{
    void mutate(boost::shared_ptr<T> ind, CmdArguments_t /*args*/) {

        int range = ind->genes_m.size();
        int position = static_cast<int>((rand() / (RAND_MAX + 1.0)) * range);
        ind->new_gene(position);
    }
};
