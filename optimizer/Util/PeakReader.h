//
// Class PeakReader
//   Implements a parser and value extractor for peak files (*.peaks)
//
// Copyright (c) 2017 - 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
// and the paper
// "Matching of turn pattern measurements for cyclotrons using multiobjective optimization"
// (https://doi.org/10.1103/PhysRevAccelBeams.22.064602)
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
#ifndef __PEAKREADER_H__
#define __PEAKREADER_H__

#include <string>
#include <map>

#include "Util/OptPilotException.h"

class PeakReader {
    
public:
    
    PeakReader(std::string filename);
    ~PeakReader();
    
    void parseFile();
    
    /**
     * @param nPeak is the peak number
     * @param radius stores result [mm]
     */
    void getPeak(int nPeak, double& radius);
    
    /**
     * @returns the number of peaks in the file
     */
    std::size_t getNumberOfPeaks();
    
private:
    /// Peak filename
    std::string filename_m;
    
    /// all found peaks < peak number, radius >
    std::map<int, double> peaks_m;
};

#endif
