//
// Namespace Options
//   The global OPAL option flags.
//   This namespace contains the global option flags.
//
// Copyright (c) 200x - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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
#ifndef OPTIONS_HH
#define OPTIONS_HH

#include "OptionTypes.h"
#include "Utilities/ClassicRandom.h"

#include <string>

namespace Options {
    /// Echo flag.
    //  If true, print an input echo.
    extern bool echo;

    /// Info flag.
    //  If true, print informative messages.
    extern bool info;
    extern int infoLevel;

    /// Trace flag.
    //  If true, print CPU time before and after each command.
    extern bool mtrace;

    /// Warn flag.
    //  If true, print warning messages.
    extern bool warn;
    extern int warnLevel;

    //  The global random generator.
    extern Random rangen;

    /// The current random seed.
    extern int seed;

    /// The frequency to dump the phase space, i.e.dump data when step%psDumpFreq==0
    extern int psDumpFreq;

    /// The frequency to dump statistical values, e.e. dump data when step%statDumpFreq==0
    extern int statDumpFreq;

    /// phase space dump flag for OPAL-cycl
    //  if true, dump phase space after each turn
    extern bool psDumpEachTurn;

    /// flag to decide in which coordinate frame the phase space will be dumped for OPAL-cycl
    //  - GLOBAL, in Cartesian frame of the global particle
    //  - BUNCH_MEAN, in Cartesian frame of the bunch mean
    //  - REFERENCE, in Cartesian frame of the reference (0) particle
    extern DumpFrame psDumpFrame;

    /// The frequency to dump single particle trajectory of particles with ID = 0 & 1
    extern int sptDumpFreq;

    /// The frequency to do particles repartition for better load balance between nodes
    extern int repartFreq;

    /// The number of bins that have to be emitted before the bin are squashed into a single bin
    extern int minBinEmitted;

    /// The number of steps into the simulation before the bins are squashed into a single bin
    extern int minStepForRebin;

    /// The frequency to reset energy bin ID for all particles
    extern int rebinFreq;

    /// The frequency to solve space charge fields.
    extern int scSolveFreq;

    // How many small timesteps are inside the large timestep used in multiple time stepping (MTS) integrator
    extern int mtsSubsteps;

    // If the distance of a particle to bunch mass larger than remotePartDel times of the rms size of the bunch in any dimension,
    // the particle will be deleted artifically to hold the accuracy of space charge calculation. The default setting of -1 stands for no deletion.
    extern double remotePartDel;

    extern bool rhoDump;

    extern bool ebDump;

    extern bool csrDump;

    // the number of refinements of the search range for the phase with maximum energy
    // if eq 0 then no autophase
    // if true opal find the phases in the cavities, such that the energy gain is at maximum
    extern int autoPhase;

    // Options for the Belos solver
    /// RCG: cycle length
    extern int numBlocks;
    /// RCG: number of recycle blocks
    extern int recycleBlocks;
    /// number of old left hand sides used to extrapolate a new start vector
    extern int nLHS;

    /// If true create symmetric distribution
    extern bool cZero;

    /// random number generator
    extern std::string rngtype;

    /// If true HDF5 files are written
    extern bool enableHDF5;

    /// If true VTK files are written
    extern bool enableVTK;

    extern bool asciidump;

    // Governs how often boundp_destroy is called to destroy lost particles
    // Mainly used in the CyclotronTracker as of now -DW
    extern int boundpDestroyFreq;

    extern double beamHaloBoundary;

    /// Do closed orbit and tune calculation only.
    extern bool cloTuneOnly;

    // Using hard edge model for calculation of path length
    extern bool idealized;

    extern bool writeBendTrajectories;

    /// opal version of input file
    extern int version;

    /// Enable AMR if true
    extern bool amr;

    /// The frequency to dump AMR grid data and particles into file
    extern int amrYtDumpFreq;
    
    /// After how many steps the AMR grid hierarchy is updated
    extern int amrRegridFreq;

    extern bool memoryDump;

    /// The constant parameter C to shift halo, by < w^4 > / < w^2 > ^2 - C (w=x,y,z)
    extern double haloShift;

    /// The frequency to delete particles (currently: OPAL-cycl only)
    extern unsigned int delPartFreq;
}

#endif // OPAL_Options_HH
