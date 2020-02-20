#ifndef OPAL_SAMPLING_METHOD_H
#define OPAL_SAMPLING_METHOD_H


#include "Sample/SampleIndividual.h"

#include <boost/smart_ptr.hpp>

#include "Comm/types.h"
#include "Util/CmdArguments.h"

class SamplingMethod
{

public:
    virtual ~SamplingMethod() {};
    virtual void create(boost::shared_ptr<SampleIndividual>& ind, size_t i) = 0;
    
    /*!
     * Allocate memory for sampling. Not every sampling method
     * requires that.
     * 
     * This function is used to reduce memory since only the
     * sampler ranks need these sampling methods.
     * 
     * @param args samler arguments
     * @param comm sampler communicator
     */
    virtual void allocate(const CmdArguments_t& /*args*/, const Comm::Bundle_t& /*comm*/) {
        /* Some sampling methods require a container.
         * In order to reduce memory only samplers should allocate
         * the memory
         */
    }
};

#endif
