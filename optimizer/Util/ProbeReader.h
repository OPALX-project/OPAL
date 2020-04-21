//
// Class ProbeReader
//   Implements a parser and value extractor for Probe loss files.
//
// Copyright (c) 2010 - 2013, Yves Ineichen, ETH Zürich
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Toward massively parallel multi-objective optimization with application to
// particle accelerators" (https://doi.org/10.3929/ethz-a-009792359)
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
#ifndef __PROBEREADER_H__
#define __PROBEREADER_H__

#include <string>
#include <vector>
#include <map>

class ProbeReader {
    
public:
    explicit ProbeReader(std::string filename);
    
    ~ProbeReader();

    void parseFile();
    
    void getVariableValue(int id, std::string varname, double& sim_value);
    
private:
    /// Probe loss filename
    std::string filename_m;
    
    /// Number of variables
    int nColumns_m;
    
    /// Number of values per variable
    int nRows_m;
    
    std::map<std::string, int> columnNamesToID_m;
    std::vector< std::vector<double> > data_m;
    
};

#endif
