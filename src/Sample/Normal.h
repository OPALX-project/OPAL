#ifndef OPAL_NORMAL_RANDOM_SAMPLING_H
#define OPAL_NORMAL_RANDOM_SAMPLING_H

#include "Sample/SamplingMethod.h"
#include "Sample/RNGStream.h"

#include <type_traits>

class Normal : public SamplingMethod
{

public:
    typedef std::normal_distribution<double> dist_t;


    Normal(double lower, double upper)
        : dist_m(0.5 * (lower + upper), (upper - lower) / 10)
        , RNGInstance_m(RNGStream::getInstance())
        , seed_m(RNGStream::getGlobalSeed())
    {}

    Normal(double lower, double upper, std::size_t seed)
        : dist_m(0.5 * (lower + upper), (upper - lower) / 10)
        , RNGInstance_m(nullptr)
        , seed_m(seed)
    {}

    ~Normal() {
        if ( RNGInstance_m)
            RNGStream::deleteInstance(RNGInstance_m);
    }

    void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) {
        ind->genes[i] = RNGInstance_m->getNext(dist_m);
    }
    
    void allocate(const CmdArguments_t& /*args*/, const Comm::Bundle_t& comm) {
        if ( !RNGInstance_m )
            RNGInstance_m = RNGStream::getInstance(seed_m + comm.island_id);
    }

private:
    dist_t dist_m;

    RNGStream *RNGInstance_m;
    
    std::size_t seed_m;
};

#endif