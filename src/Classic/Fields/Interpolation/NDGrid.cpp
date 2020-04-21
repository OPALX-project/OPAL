/* 
 *  Copyright (c) 2017, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met: 
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer. 
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to 
 *     endorse or promote products derived from this software without specific 
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <math.h>
#include <iomanip>
#include "Utilities/OpalException.h"
#include "Fields/Interpolation/Mesh.h"
#include "Fields/Interpolation/NDGrid.h"

namespace interpolation {
/// ////// NDGrid ///////
NDGrid::NDGrid() : coord_m(), maps_m(), constantSpacing_m(false)  {
}

NDGrid::NDGrid(std::vector<int> size, std::vector<const double *> gridCoordinates)
    : coord_m(), maps_m(), constantSpacing_m(false) { 
    for (unsigned int i=0; i<size.size(); i++) {
        if (size[i] < 1) {
            throw(OpalException("NDGrid::NDGrid(...)",
                                "ND Grid must be at least 1x1x...x1 grid"));
        }
        coord_m.push_back(std::vector<double>(gridCoordinates[i],
                                              gridCoordinates[i] + size[i])); 
    }
    setConstantSpacing();
}

NDGrid::NDGrid(int nDimensions, int* size, double* spacing, double* min)
    : coord_m(nDimensions), maps_m(), constantSpacing_m(true) {
    for (int i=0; i<nDimensions; i++) {
        if (size[i] < 1) {
            throw(OpalException("NDGrid::NDGrid(...)",
                                "ND Grid must be at least 1x1x...x1 grid"));
        }
        coord_m[i] = std::vector<double>(size[i]);
        for (unsigned int j=0; j<coord_m[i].size(); j++) {
            coord_m[i][j] = min[i] + j*spacing[i];
        }
    }
}

NDGrid::NDGrid(std::vector< std::vector<double> > gridCoordinates)
    : coord_m(gridCoordinates), maps_m(), constantSpacing_m(false) {
    for (unsigned int i=0; i<gridCoordinates.size(); i++) {
        if (gridCoordinates[i].size() < 1) {
            throw (OpalException("NDGrid::NDGrid(...)",
                                "ND Grid must be at least 1x1x...x1 grid"));
        }
    }
    setConstantSpacing();
}


double* NDGrid::newCoordArray  ( const int& dimension)  const {
    double * array = new double[coord_m[dimension].size() ];
    for (unsigned int i=0; i<coord_m[dimension].size(); i++) {
        array[i] = coord_m[dimension][i];
    }
    return array;
}

//Mesh::Iterator wraps around a std::vector<int>
//it[0] is least significant, it[max] is most signifcant
Mesh::Iterator& NDGrid::addEquals(Mesh::Iterator& lhs, int difference) const {
    if (difference < 0) {
        return subEquals(lhs, -difference);
    }
    std::vector<int> index(coord_m.size(),0);
    std::vector<int> content(coord_m.size(),1);
    for (int i = int(index.size()-2); i >= 0; i--) {
        content[i] = content[i+1]*coord_m[i+1].size(); //content could be member variable
    }
    for (int i = 0; i < int(index.size()); i++) {
        index[i] = difference/content[i];
        difference -= index[i] * content[i];
    }
    for (unsigned int i=0; i<index.size(); i++) {
        lhs.state_m[i] += index[i];
    }
    for (int i = int(index.size())-1; i > 0; i--) {
        if (lhs.state_m[i] > int(coord_m[i].size())) {
            lhs.state_m[i-1]++;
            lhs.state_m[i] -= coord_m[i].size();
        }
    }

    return lhs;
}

Mesh::Iterator& NDGrid::subEquals(Mesh::Iterator& lhs, int difference) const {
    if (difference < 0) {
        return addEquals(lhs, -difference);
    }
    std::vector<int> index(coord_m.size(),0);
    std::vector<int> content(coord_m.size(),1);
    for (int i = int(index.size()-2); i >= 0; i--) {
        content[i] = content[i+1]*coord_m[i+1].size(); //content could be member variable
    }
    for (int i = 0; i < int(index.size()); i++) {
        index[i]    = difference/content[i];
        difference -= index[i] * content[i];
    }
    for (unsigned int i = 0; i < index.size(); i++) {
        lhs.state_m[i] -= index[i];
    }
    for (int i=int(index.size())-1; i>0; i--) {
        if (lhs.state_m[i] < 1) {
            lhs.state_m[i-1]--;
            lhs.state_m[i] += coord_m[i].size();
        }
    }
    return lhs;
}

Mesh::Iterator& NDGrid::addEquals(Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const {
    return addEquals(lhs, rhs.toInteger());
}

Mesh::Iterator& NDGrid::subEquals(Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const {
    return subEquals(lhs, rhs.toInteger());
}

Mesh::Iterator& NDGrid::addOne(Mesh::Iterator& lhs) const {
    int i = coord_m.size()-1;
    while (lhs[i] == int(coord_m[i].size()) && i>0) {
        lhs[i]=1; i--;
    }
    lhs[i]++;
    return lhs;
}

Mesh::Iterator& NDGrid::subOne(Mesh::Iterator& lhs) const {
    lhs[coord_m.size()-1] -= 1;

    int i = coord_m.size()-1;
    while (lhs[i] == 0 && i>0) {
        lhs.state_m[i] = coord_m[i].size();
        i--;
        lhs.state_m[i]--;
    }
    return lhs;
}

void NDGrid::setConstantSpacing(double tolerance_m) {
    constantSpacing_m = true;
    for (unsigned int i = 0; i < coord_m.size(); i++) {
        for (unsigned int j = 0; j < coord_m[i].size()-1; j++) {
            double coord_j1 = coord_m[i][j+1];
            double coord_j0 = coord_m[i][j];
            double coord_1 = coord_m[i][1];
            double coord_0 = coord_m[i][0];
            if (std::abs(1-(coord_j1-coord_j0)/(coord_1-coord_0)) > tolerance_m) {
                constantSpacing_m = false;
                return;
            }
        }
    }
}

bool NDGrid::isGreater(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const {
    unsigned int i = 0;
    // if all equal; rhs[i] = rhs.last
    while (i < rhs.state_m.size()-1 && lhs.state_m[i] == rhs.state_m[i]) {
        i++; 
    }
    return (lhs[i] > rhs[i]);
}

int NDGrid::toInteger(const Mesh::Iterator& lhs) const {
    int difference = 0;
    std::vector<int> index  (coord_m.size(),0);
    std::vector<int> content(coord_m.size(),1);
    for (int i = int(index.size()-2); i >= 0; i--) {
        content[i]  = content[i+1]*coord_m[i+1].size();
    }
    for (int i = 0; i < int(index.size()); i++) {
        difference += (lhs.state_m[i]-1) * (content[i]);
    }
    return difference;
}

Mesh::Iterator NDGrid::getNearest(const double* position) const {
    std::vector<int>    index(coord_m.size());
    std::vector<double> pos(position, position+coord_m.size());
    lowerBound(pos, index);
    for (unsigned int i = 0; i < coord_m.size(); i++)
    {
        if (index[i] < int(coord_m[i].size()-1) && index[i] >= 0) {
            index[i] += (2*(position[i] - coord_m[i][index[i]])  > coord_m[i][index[i]+1]-coord_m[i][index[i]] ? 2 : 1);
        } else {
            index[i]++;
        }
        if (index[i] < 1) {
            index[i] = 1; 
        }
        if (index[i] > int(coord_m[i].size())) {
            index[i] = coord_m[i].size(); 
        }
    }
    return Mesh::Iterator(index, this);
}


Mesh* NDGrid::dual() const {
    std::vector<std::vector<double> > coord(coord_m.size());
    for (size_t i = 0; i < coord.size(); ++i) {
        if (coord_m[i].size() <= 1) {
            throw(OpalException("NDGrid::dual(...)",
                                "ND Grid must be at least 2x2x...x2 grid"));
        }
        coord[i] = std::vector<double>(coord_m[i].size()-1);
        for (size_t j = 0; j < coord[i].size(); ++j) {
            coord[i][j] = (coord_m[i][j] + coord_m[i][j+1])/2.;
        }
    }
    return new NDGrid(coord);
}


////// NDGrid END ////////
}
