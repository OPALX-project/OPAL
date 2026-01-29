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
#include "ProbeReader.h"

#include "Util/OptPilotException.h"

#include <cstring>
#include <fstream>

ProbeReader::ProbeReader(std::string filename) :
    filename_m(filename),
    nColumns_m(0),
    nRows_m(0),
    data_m(0)
{ }


ProbeReader::~ProbeReader() { }


void ProbeReader::parseFile() {
    
    nColumns_m = 0;
    nRows_m = 0;
    
    std::ifstream probe;
    
    probe.open(filename_m.c_str(), std::ios::in);
    if ( !probe ) {
        throw OptPilotException("ProbeReader::parseFile()",
                                "Error opening file " + filename_m);
    }
    
    std::string header;
    std::getline(probe, header, '\n');
    
    if( header.find("# Element\n") == std::string::npos ) {
        throw OptPilotException("ProbeReader::parseFile()",
                                "Error parsing Probe header!");
    }
    
    char *token = std::strtok(&header[0], " ");
    
    // parse header
    while ( token != NULL ) {
        
        if ( std::string(token) == "Element")
            token = std::strtok(NULL, " ");  // skip name
        else if ( token[0] != ')' && token[0] != '(' && token[0] != '#') {
            std::string varname = std::string(token);
            
            if ( varname.back() == ',' )
                varname.pop_back();
            
            columnNamesToID_m[varname] = nColumns_m++;
        }
        token = std::strtok(NULL, " ");
    }
    
    // parse values
    data_m.resize(nColumns_m);
    
    std::string line;
    while ( std::getline(probe, line) ) {
        
        ++nRows_m;
        
        token = std::strtok(&line[0], " ");
        
        // skip first (probe name)
        token = std::strtok(NULL, " ");
        
        int i = 0;
        
        while ( token != NULL ) {
            data_m[i++].push_back( std::atof(token) );
            token = std::strtok(NULL, " ");
        }
    }
    
    probe.close();
}


void ProbeReader::getVariableValue(int id, std::string varname, double& sim_value) {
    
    
    int varindex = 0;
    if(columnNamesToID_m.count(varname) > 0) {
        varindex = columnNamesToID_m[varname];
    } else {
        throw OptPilotException("ProbeReader::getVariableValue",
                                "variable name!");
    }
    
    int col = 0;
    if(columnNamesToID_m.count("id") > 0) {
        col = columnNamesToID_m["id"];
    } else {
        throw OptPilotException("ProbeReader::getVariableValue",
                                "ID variable not found!");
    }
    
    int row = -1;
    for (unsigned int i = 0; i < data_m[col].size(); ++i) {
        if ( data_m[col][i] == id ) {
            row = i;
            break;
        }
    }
    
    if ( row < 0 )
        throw OptPilotException("ProbeReader::getVariableValue",
                                "Appropriate value not found!");
    
    sim_value = data_m[varindex][row];
}
