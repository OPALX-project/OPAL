#ifndef __OPAL_SAMPLER_H__
#define __OPAL_SAMPLER_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <fstream>


#include "Comm/types.h"
#include "Util/Types.h"
#include "Util/CmdArguments.h"

#include "Optimizer/Optimizer.h"
#include "Sample/SampleIndividual.h"
#include "Sample/SamplingMethod.h"

#include <boost/smart_ptr.hpp>
#include <boost/chrono.hpp>

#include <memory>

#include <queue>

#include <boost/property_tree/ptree.hpp>


/**
 *  \class Sampler
 *  \brief Implementing sampling
 */
// template<
//       template <class> class SamplingOperator
// >
class Sampler : public Optimizer {

public:

    /** This constructor should never be called.
     *  It's provided due to the inheritance of SamplePilot from
     *  Pilot
     *
     */
    Sampler(Expressions::Named_t objectives,
            Expressions::Named_t constraints,
            DVarContainer_t dvars,
            size_t dim, Comm::Bundle_t comms,
            CmdArguments_t args,
            std::vector<double> hypervolRef,
            int nrWorkerGroups);


    /**
     *  Retrieves all (for the sampler) relevant arguments specified on the
     *  command line, initializes the variator and sets up statistics and
     *  debug traces.
     *
     *  @param[in] sampleMethods per design variable (dvar)
     *  @param[in] dvars of sampling
     *  @param[in] comms available to the sampler
     *  @param[in] args the user passed on the command line
     */
    Sampler(const std::map< std::string,
                            std::shared_ptr<SamplingMethod>
                >& sampleMethods,
            Expressions::Named_t objectives,
            DVarContainer_t dvars,
            Comm::Bundle_t comms,
            CmdArguments_t args);

    /// Initialization and start algorithm
    virtual void initialize();

    /// type used in solution state exchange with other optimizers
    typedef std::vector< SampleIndividual > SolutionState_t;

protected:

    // implementing poller hooks
    bool onMessage(MPI_Status status, size_t length);
    void postPoll();

    void setupPoll() {}
    void prePoll() {}
    void onStop() {}

    // helper sending evaluation requests to the pilot
    void dispatch_forward_solves();

private:

    std::map<std::string,
             std::shared_ptr<SamplingMethod>
        > sampleMethods_m;

    // global index (for job id)
    int gid;

    int my_local_pid_;

    typedef SampleIndividual  Individual_t;

    /// communicator bundle for the optimizer
    Comm::Bundle_t comms_;

    /// mapping from unique job ID to individual
    std::map<size_t, boost::shared_ptr<Individual_t> > jobmapping_m;

    std::queue<boost::shared_ptr<Individual_t> > individuals_m;

    /// bounds on each specified gene
    bounds_t dVarBounds_m;

    /// design variables
    DVarContainer_t dvars_m;

    /// objectives
    Expressions::Named_t objectives_m;

    int nSamples_m;


    /// command line arguments specified by the user
    CmdArguments_t args_;

    /// current generation
    int act_sample_m;

    int done_sample_m;

    enum State {
        SUBMIT,
        STOP,
        TERMINATE
    };

    State curState_m;

    /// Dumps id, design variables and bound
    std::string resultFile_m;
    std::string resultDir_m;

    boost::property_tree::ptree tree_m;

    void dumpIndividualsToJSON_m();

    void addIndividualToJSON_m(const boost::shared_ptr<Individual_t>& ind);

    void runStateMachine();

    void createNewIndividual_m();
};

#endif