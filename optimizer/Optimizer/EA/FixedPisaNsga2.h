#ifndef __FIXED_PISA_NSGA2_H__
#define __FIXED_PISA_NSGA2_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <deque>
#include <fstream>
#include <map>
#include <string>
#include <sstream>
#include <utility>
#include <vector>

#include "Comm/types.h"
#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Util/Statistics.h"

#include "Optimizer/Optimizer.h"
#include "Optimizer/EA/Individual.h"
#include "Optimizer/EA/Variator.h"

#include <boost/smart_ptr.hpp>
#include <boost/chrono.hpp>
#include <boost/property_tree/ptree.hpp>

#include "Util/Trace/Trace.h"


/**
 *  \class FixedPisaNsga2
 *  \brief Implementing the Variator for the PISA state machine.
 *
 *  @see http://www.tik.ee.ethz.ch/pisa/
 *
 *  The convergence behavior of the optimizer can be steered in 3 ways,
 *  corresponding command line arguments are given in brackets:
 *    - limit the number of generations (maxGenerations),
 *    - specify a target hypervolume (expected-hypervol) and tolerance
 *      (epsilon)
 *    - specify a minimal hypervolume progress (conv-hvol-prog), relative to
 *      the last generation, ((prev - new)/prev) that has to be attained to
 *      continue optimizing.
 */
template<
      template <class> class CrossoverOperator
    , template <class> class MutationOperator
>
class FixedPisaNsga2 : public Optimizer {

public:

    /**
     *  Retrieves all (for the optimizer) relevant arguments specified on the
     *  command line, initializes the variator and sets up statistics and
     *  debug traces.
     *
     *  @param[in] objectives of optimization problem
     *  @param[in] constraints of optimization problem
     *  @param[in] dvars of optimization problem
     *  @param[in] dim number of objectives
     *  @param[in] comms available to the optimizer
     *  @param[in] args the user passed on the command line
     *  @param[in] hypervolRef hypervolume reference point
     */
    FixedPisaNsga2(Expressions::Named_t objectives,
                   Expressions::Named_t constraints,
                   DVarContainer_t dvars, size_t dim, Comm::Bundle_t comms,
                   CmdArguments_t args,
                   std::vector<double> hypervolRef,
                   int nrWorkerGroups);

    ~FixedPisaNsga2();

    /// Starting selection algorithm and variator PISA state machine
    virtual void initialize();

    /// type used in solution state exchange with other optimizers
    typedef std::vector< Individual > SolutionState_t;
    /// type of our variator
    typedef Individual Individual_t;
    typedef Variator< Individual_t, CrossoverOperator, MutationOperator >
        Variator_t;
    /// alias for usage in template
    using individual = typename FixedPisaNsga2::Individual_t;


protected:

    /// Write the variator config file
    void writeVariatorCfg();

    ///@{ implementing poller hooks
    virtual bool onMessage(MPI_Status status, size_t length);
    virtual void postPoll();

    virtual void setupPoll() {}
    virtual void prePoll() {
        // std::ostringstream debug;
        // debug << "IN PRE POLL: ";
        // debug << getStateString(curState_m) << std::endl;
        // progress_->log(debug);
    }
    virtual void onStop() {}
    ///@}
    /// helper sending evaluation requests to the pilot
    void dispatch_forward_solves();


private:

    int my_local_pid_;

    /// all PISA states
    enum PisaState_t {
          Initialize         = 0
        , InitializeSelector = 1
        , Variate            = 2
        , Select             = 3
        , Stop               = 4
        , VariatorStopped    = 5
        , VariatorTerminate  = 6
        /* , SelectorStopped    = 7 */
        /* , Reset              = 8 */
        /* , ReadyForReset      = 9 */
        /* , ReadyForResetS     = 10 */
        /* , Restart            = 11 */
    };

    std::string getStateString(PisaState_t) const;

    // selector parameters
    int seed;   /* seed for random number generator */
    const int tournament_m = 1;  /* number of opponents for mating selection */

    /// the current state of the state machine
    PisaState_t curState_m;

    /// collect some statistics of rejected and accepted individuals
    boost::scoped_ptr<Statistics<size_t> > statistics_;

    boost::scoped_ptr<Variator_t> variator_m;

    std::vector<unsigned int> pp_all;        ///< IDs of population
    std::vector<unsigned int> parent_queue_; ///< IDs that will make new offspring
    // std::set<unsigned int> archive_;
    std::set<unsigned int> to_selector_;     ///< Successfully run IDs to go into population

    // to compute the front
    std::vector<int> copies; ///< number of individuals in the n-th front
    std::vector<double> dist;
    std::vector< std::vector<int> > front; ///< individuals in each front
    std::map<size_t, double> fitness_; ///< map between id and fitness (sum of front number and dist)

    /// communicator bundle for the optimizer
    Comm::Bundle_t comms_;

    /// buffer holding all finished job id's
    std::deque<unsigned int> finishedBuffer_m;

    /// mapping from unique job ID to individual
    std::map<size_t, boost::shared_ptr<Individual> > jobmapping_m;

    /// indicating if initial population has been created
    bool initialized_m;

    /// bounds on each specified gene
    bounds_t dVarBounds_m;
    /// objectives
    Expressions::Named_t objectives_m;
    /// constraints
    Expressions::Named_t constraints_m;
    /// design variables
    DVarContainer_t dvars_m;

    /// command line arguments specified by the user
    CmdArguments_t args_m;

    /// size of initial population
    size_t alpha_m;
    /// initial population optimization flag (increases initial population)
    bool initialOptimization_m;
    /// enforce strict population size
    bool birthControl_m;
    /// population file to be started from
    std::string file_start_m;

    /// number of parents the selector chooses
    //size_t mu_m;
    /// number of children the variator produces
    size_t lambda_m;
    /// number of objectives
    size_t dim_m;
    /// current generation
    size_t act_gen = 1;
    /// maximal generation (stopping criterion)
    size_t maxGenerations_m;
    /// number of individuals running
    int num_workergroups_m;

    /// result file name
    std::string resultFile_m;
    std::string resultDir_m;


    // dump frequency
    int dump_freq_m;
    /// dump offspring / parents flag
    bool dump_offspring_m;
    /// convergence accuracy if maxGenerations not set
    double hvol_eps_;
    double expected_hvol_;
    double current_hvol_;
    double conv_hvol_progress_;
    double hvol_progress_;

    /// hypervolume reference point
    std::vector<double> hvol_ref_m;

    /// file header for result files contains this parameter description
    std::string file_param_descr_;

    boost::chrono::system_clock::time_point run_clock_start_;
    boost::chrono::system_clock::time_point last_clock_;

    // DEBUG output helpers
    boost::scoped_ptr<Trace> job_trace_;
    boost::scoped_ptr<Trace> progress_;


    // entry point for starting the selector side of the PISA state machine
    void startSelector(std::string filename_base);

    /// executes one loop of the PISA state machine
    void runStateMachine();

    /// passes finished individuals to the selector
    void toSelectorAndCommit();

    /// how often do we exchange solutions with other optimizers
    size_t exchangeSolStateFreq_m;

    /// if necessary exchange solution state with other optimizers
    void exchangeSolutionStates();

    // Selector methods
    void selection();
    void mergeOffspring();
    void calcFitnesses();
    void calcDistances();
    void environmentalSelection();
    void matingSelection();
    int dominates(unsigned int p_ind_a, unsigned int p_ind_b);

    /// Dumps index, objective values and bit string of all individuals in
    /// global_population.
    void dumpPopulationToFile();
    void dumpPopulationToJSON();
    void dumpIndividualToFile(int id,
                              boost::shared_ptr<individual>& ind,
                              std::ofstream& file,
                              const size_t numDigits);
    void dumpIndividualToJSON(int id,
                              boost::shared_ptr<individual>& ind,
                              boost::property_tree::ptree& tree);

    /**
     *  Get a random integer between [0, range]
     *  @param[in] range of random number
     *  @return random integer value between [0, range]
     */
    int irand(int range) {
        return (int) ((double) range * (double) rand() / (RAND_MAX + 1.0));
    }
};

#include "Optimizer/EA/FixedPisaNsga2.tcc"

#endif