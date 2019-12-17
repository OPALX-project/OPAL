#include "OptionTypes.h"
#include "Utilities/ClassicRandom.h"
#include <string>

namespace Options {
    // The global program options.
    bool echo = false;
    bool info = true;
    int infoLevel = 1;
    bool warn = true;
    int warnLevel = 1;
    bool csrDump = false;
    bool ppdebug = false;

    // If true create symmetric distribution
    bool cZero = false;

    bool asciidump = false;

    // If the distance of a particle to bunch mass larger than remotePartDel times of the rms size of the bunch in any dimension,
    // the particle will be deleted artifically to hold the accuracy of space charge calculation. The default setting of -1 stands for no deletion.
    double remotePartDel = 0.0;

    double beamHaloBoundary = 0;

    bool writeBendTrajectories = false;

    // The global program options.
    bool mtrace = false;
    bool psDumpEachTurn = false;
    DumpFrame psDumpFrame = GLOBAL;
    bool rhoDump = false;
    bool ebDump = false;

    bool enableHDF5 = true;

    // The global random generator.
    Random rangen;

    // The current random seed.
    int seed  = 123456789;

    // the number of refinements of the search range for the phase with maximum energy
    // if eq 0 then no autophase
    int autoPhase = 6;

    // The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0
    int psDumpFreq = 10;
    // // The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0
    // double rDump = 0.0;

    // The frequency to dump statistical quantities such as beam RMS properties, i.e. dump
    // when step%statDumpFreq == 0.
    int statDumpFreq = 10;

    // The frequency to dump single particle trajectory of particles with ID = 0 & 1
    int sptDumpFreq = 1;

    // The frequency to do particles repartition for better load balance between nodes
    int repartFreq = 10;

    // The number of bins that have to be emitted before the bin are squashed into a single bin
    int minBinEmitted = 10;

    // The number of steps into the simulation before the bins are squashed into a single bin
    int minStepForRebin = 200;

    // The frequency to reset energy bin ID for all particles
    int rebinFreq = 100;

    /// The frequency to solve space charge fields.
    int scSolveFreq = 1;

    // How many small timesteps are inside the large timestep used in multiple time stepping (MTS) integrator
    int mtsSubsteps = 1;

    // The frequency to dump the particle-geometry surface interation data, -1 stands for no dump.
    int surfDumpFreq = -1;

    // Options for the Belos solver
    int numBlocks = 0;
    int recycleBlocks = 0;
    int nLHS = 1;

    std::string rngtype = std::string("RANDOM");

    bool cloTuneOnly = false;

    // Governs how often boundp_destroy is called to destroy lost particles
    // Mainly used in the CyclotronTracker as of now -DW
    int boundpDestroyFreq = 10;

    // Using hard edge model for calculation of path length
    bool idealized = false;

    // opal version of input file
    int version = 10000;

#ifdef ENABLE_AMR
    bool amr = false;

    /// The frequency to dump AMR grid data and particles into file
    int amrYtDumpFreq = 10;
    
    int amrRegridFreq = 10;
#endif

    bool memoryDump = false;

    double haloShift = 0.0;

    unsigned int delPartFreq = 1;
}