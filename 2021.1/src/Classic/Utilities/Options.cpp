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
#include "OptionTypes.h"
#include "Utilities/ClassicRandom.h"

#include <string>

namespace Options {
    // The global program options.
    bool echo = false;

    bool info = true;
    int infoLevel = 1;

    bool mtrace = false;

    bool warn = true;
    int warnLevel = 1;

    Random rangen;
    int seed  = 123456789;

    int psDumpFreq = 10;

    int statDumpFreq = 10;

    bool psDumpEachTurn = false;

    DumpFrame psDumpFrame = GLOBAL;

    int sptDumpFreq = 1;

    int repartFreq = 10;

    int minBinEmitted = 10;

    int minStepForRebin = 200;

    int rebinFreq = 100;

    int scSolveFreq = 1;

    int mtsSubsteps = 1;

    double remotePartDel = 0.0;

    bool rhoDump = false;

    bool ebDump = false;

    bool csrDump = false;

    int autoPhase = 6;

    int numBlocks = 0;

    int recycleBlocks = 0;

    int nLHS = 1;

    bool cZero = false;

    std::string rngtype = std::string("RANDOM");

    bool enableHDF5 = true;

    bool enableVTK = true;

    bool asciidump = false;

    int boundpDestroyFreq = 10;

    double beamHaloBoundary = 0;

    bool cloTuneOnly = false;

    bool idealized = false;

    bool writeBendTrajectories = false;

    int version = 10000;

    bool amr = false;

    int amrYtDumpFreq = 10;

    int amrRegridFreq = 10;

    bool memoryDump = false;

    double haloShift = 0.0;

    unsigned int delPartFreq = 1;
}
