#ifndef __OPALSIMULATION_H__
#define __OPALSIMULATION_H__

#include <string>
#include <map>

#include "mpi.h"

#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Simulation/Simulation.h"
#include "Simulation/GenerateOpalSimulation.h"

#include "boost/smart_ptr.hpp"

/**
 *  \class OpalSimulation
 *  \brief Concrete implementation of an Opal simulation wrapper.
 *
 *  @see http://amas.web.psi.ch
 *
 *  In order to work properly the user must ensure that the environment
 *  variables
 *
 *    - SIMTMPDIR (temporary directory for simulation input data) and
 *    - TEMPLATES (directory containing tmpl file)
 *
 *  are specified correctly.
 */
class OpalSimulation : public Simulation {

public:

    /**
     *  Setup OPAL run.
     *
     *  @param[in] objectives of the optimization problem
     *  @param[in] constraints of the optimization problem
     *  @param[in] params
     *  @param[in] name of the simulation
     *  @param[in] comm MPI communicator used for running the simulation
     *  @param[in] args command line arguments passed to the framework
     */
    OpalSimulation(Expressions::Named_t objectives,
                   Expressions::Named_t constraints,
                   Param_t params, std::string name, MPI_Comm comm,
                   CmdArguments_t args);

    virtual ~OpalSimulation();

    /// Calls Opal through Opal-lib wrapper and returns when simulation has
    /// either failed or finished.
    void run();

    /// Parse SDDS stat file and build up requested variable dictionary.
    void collectResults();

    /// returns container containing all requested variables with results
    reqVarContainer_t getResults() { return requestedVars_; }
    
    /// set job id (SAMPLE command)
    void setFilename(int id) { id_m = id; }

    std::string getSimDirectory() { return simDir_; }

    /// remove temporary simulation files (if Boost filesystem library is not
    /// available, this call do nothing).
    void cleanUp();

    void cleanUp(const std::vector<std::string>& keep);

private:

    /// identification of the simulation (corresponding to output filename)
    std::string simulationName_;
    /// full path of simulation directory (where simulation will be run)
    std::string simulationDirName_;
    /// temporary directory for simulation data (environment var SIMTMPDIR)
    std::string simTmpDir_;
    std::string simDir_;

    /// holds current directory (for restoring)
    std::string pwd_;

    /// stream buffer to redirect output
    std::streambuf* strm_buffer_;
    /// stream buffer to redirect stderr
    std::streambuf* strm_err_;

    /// variable dictionary holding requested optimizer values
    std::map<std::string, std::string> userVariables_;

    /// holds solutions returned to the optimizer
    reqVarContainer_t requestedVars_;

    Expressions::Named_t objectives_;
    Expressions::Named_t constraints_;

    MPI_Comm comm_;

    /// object to generate simulation input files
    boost::scoped_ptr<GenerateOpalSimulation> gs_;

    /// job id (SAMPLE command)
    int id_m;

    /// mark a solution as invalid
    void invalidBunch();

    /// check if we already have simulated the current set of design vars
    bool hasResultsAvailable();
    /// create directories, input files, fieldmaps...
    void setupSimulation();

    /// get variables for expression evaluation from SDDS file. If failed returns false
    bool getVariableDictionary(variableDictionary_t& dictionary,
                               const std::string& filename,
                               const Expressions::Expr_t* const expression);

    /// redirect stdout and stderr to file
    void redirectOutToFile();
    /// restore stdout and stderr to default
    void restoreOut();
};

#endif