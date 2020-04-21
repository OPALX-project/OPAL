/*
 *  Copyright (c) 2012, Chris Rogers
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

#include "Fields/Interpolation/ThreeDGrid.h"
#include "Utilities/LogicalError.h"

namespace interpolation {

ThreeDGrid::ThreeDGrid()
    : x_m(2, 0), y_m(2, 0), z_m(2, 0), xSize_m(0), ySize_m(0), zSize_m(0),
      maps_m(0), constantSpacing_m(false) {
    setConstantSpacing();
    x_m[1] = y_m[1] = z_m[1] = 1.;
}

ThreeDGrid::ThreeDGrid(int xSize, const double *x,
                       int ySize, const double *y,
                       int zSize, const double *z)
    : x_m(x, x+xSize), y_m(y, y+ySize), z_m(z, z+zSize),
      xSize_m(xSize), ySize_m(ySize), zSize_m(zSize),
      maps_m(), constantSpacing_m(false) {
  setConstantSpacing();
  if (xSize_m < 2 || ySize_m < 2 || zSize_m < 2)
    throw(LogicalError(
        "ThreeDGrid::ThreeDGrid(...)",
        "3D Grid must be at least 2x2x2 grid"
        )
    );
}

ThreeDGrid::ThreeDGrid(std::vector<double> x,
                       std::vector<double> y,
                       std::vector<double> z)
    : x_m(x), y_m(y), z_m(z),
      xSize_m(x.size()), ySize_m(y.size()), zSize_m(z.size()),
      maps_m(), constantSpacing_m(false) {
  setConstantSpacing();
  if (xSize_m < 2 || ySize_m < 2 || zSize_m < 2)
      throw(LogicalError(
          "ThreeDGrid::ThreeDGrid(...)",
          "3D Grid must be at least 2x2x2 grid"
          )
      );
}

ThreeDGrid::ThreeDGrid(double dX, double dY, double dZ,
                       double minX, double minY, double minZ,
                       int numberOfXCoords, int numberOfYCoords,
                       int numberOfZCoords)
    :  x_m(numberOfXCoords), y_m(numberOfYCoords), z_m(numberOfZCoords),
       xSize_m(numberOfXCoords), ySize_m(numberOfYCoords),
       zSize_m(numberOfZCoords), maps_m(), constantSpacing_m(true) {
    for (int i = 0; i < numberOfXCoords; i++)
        x_m[i] = minX+i*dX;
    for (int j = 0; j < numberOfYCoords; j++)
        y_m[j] = minY+j*dY;
    for (int k = 0; k < numberOfZCoords; k++)
        z_m[k] = minZ+k*dZ;

    setConstantSpacing();
}

ThreeDGrid::~ThreeDGrid() {
}

// state starts at 1,1,1
Mesh::Iterator ThreeDGrid::begin() const {
    return Mesh::Iterator(std::vector<int>(3, 1), this);
}

Mesh::Iterator ThreeDGrid::end() const {
    std::vector<int> end(3, 1);
    end[0] = xSize_m+1;
    return Mesh::Iterator(end, this);
}

void ThreeDGrid::getPosition(const Mesh::Iterator& it,
                             double * position) const {
    position[0] = x(it.state_m[0]);
    position[1] = y(it.state_m[1]);
    position[2] = z(it.state_m[2]);
}

Mesh* ThreeDGrid::dual() const {
    std::vector<double> new_x(x_m.size()-1);
    std::vector<double> new_y(y_m.size()-1);
    std::vector<double> new_z(z_m.size()-1);
    for (size_t i = 0; i < x_m.size()-1; ++i) {
        new_x[i] = (x_m[i]+x_m[i+1])/2.;
    }
    for (size_t i = 0; i < y_m.size()-1; ++i) {
        new_y[i] = (y_m[i]+y_m[i+1])/2.;
    }
    for (size_t i = 0; i < z_m.size()-1; ++i) {
        new_z[i] = (z_m[i]+z_m[i+1])/2.;
    }
    return new ThreeDGrid(new_x, new_y, new_z);
}

Mesh::Iterator& ThreeDGrid::addEquals
    (Mesh::Iterator& lhs, int difference) const {
    if (difference < 0)
        return subEquals(lhs, -difference);

    lhs.state_m[0] +=  difference/(ySize_m*zSize_m);
    difference     = difference%(ySize_m*zSize_m);
    lhs.state_m[1] += difference/(zSize_m);
    lhs.state_m[2] += difference%(zSize_m);

    if (lhs.state_m[1] > ySize_m) {
        lhs.state_m[0]++;
        lhs.state_m[1] -= ySize_m;
    }
    if (lhs.state_m[2] > zSize_m) {
        lhs.state_m[1]++;
        lhs.state_m[2] -= zSize_m;
    }

    return lhs;
}

Mesh::Iterator& ThreeDGrid::subEquals
                                  (Mesh::Iterator& lhs, int difference) const {
    if (difference < 0)
        return addEquals(lhs, -difference);

    lhs.state_m[0] -= difference/(ySize_m*zSize_m);
    difference     = difference%(ySize_m*zSize_m);
    lhs.state_m[1] -= difference/(zSize_m);
    lhs.state_m[2] -= difference%(zSize_m);

    while (lhs.state_m[2] < 1) {
        lhs.state_m[1]--;
        lhs.state_m[2] += zSize_m;
    }
    while (lhs.state_m[1] < 1) {
        lhs.state_m[0]--;
        lhs.state_m[1] += ySize_m;
    }

    return lhs;
}

void ThreeDGrid::vectorLowerBound(std::vector<double> vec,
                                  double x,
                                  int& index) {
    if (x < vec[0]) {
        index = -1;
        return;
    }
    if (x >= vec.back()) {
        index = vec.size()-1;
        return;
    }
    size_t xLower = 0;
    size_t xUpper = vec.size()-1;
    while (xUpper - xLower > 1) {
        index = (xUpper+xLower)/2;
        if (x >= vec[index]) {
            xLower = index;
        } else {
            xUpper = index;
        }
    }
    index = xLower;
}

Mesh::Iterator& ThreeDGrid::addEquals
                      (Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const {
    return addEquals(lhs, rhs.toInteger());
}

Mesh::Iterator& ThreeDGrid::subEquals
                        (Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const {
    return subEquals(lhs, rhs.toInteger());
}

Mesh::Iterator& ThreeDGrid::addOne(Mesh::Iterator& lhs) const {
    if (lhs.state_m[1] == ySize_m && lhs.state_m[2] == zSize_m) {
        ++lhs.state_m[0];
        lhs.state_m[1] = lhs.state_m[2] = 1;
    } else if (lhs.state_m[2] == zSize_m) {
        ++lhs.state_m[1];
        lhs.state_m[2] = 1;
    } else {
        ++lhs.state_m[2];
    }
    return lhs;
}

Mesh::Iterator& ThreeDGrid::subOne(Mesh::Iterator& lhs) const {
    if (lhs.state_m[1] == 1 && lhs.state_m[2] == 1) {
        --lhs.state_m[0];
        lhs.state_m[1] = ySize_m;
        lhs.state_m[2] = zSize_m;
    } else if (lhs.state_m[2] == 1) {
        --lhs.state_m[1];
        lhs.state_m[2] = zSize_m;
    } else {
        --lhs.state_m[2];
    }
    return lhs;
}

// Check relative position
bool ThreeDGrid::isGreater
                (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const {
    if (lhs.state_m[0] > rhs.state_m[0])
        return true;
    else if (lhs.state_m[0] == rhs.state_m[0] &&
             lhs.state_m[1] > rhs.state_m[1])
        return true;
    else if (lhs.state_m[0] == rhs.state_m[0] &&
             lhs.state_m[1] == rhs.state_m[1] &&
             lhs.state_m[2] > rhs.state_m[2])
        return true;
    return false;
}

// remove *map if it exists; delete this if there are no more VectorMaps
void ThreeDGrid::remove(VectorMap* map) {
    std::vector<VectorMap*>::iterator it =
                                   std::find(maps_m.begin(), maps_m.end(), map);
    if (it < maps_m.end()) {
        maps_m.erase(it);
    }
    if (maps_m.begin() == maps_m.end()) {
        delete this;
    }
}

// add *map if it has not already been added
void ThreeDGrid::add(VectorMap* map) {
    std::vector<VectorMap*>::iterator it =
                                   std::find(maps_m.begin(), maps_m.end(), map);
    if (it == maps_m.end()) {
        maps_m.push_back(map);
    }
}

void ThreeDGrid::setConstantSpacing() {
    constantSpacing_m = true;
    for (unsigned int i = 0; i < x_m.size()-1; i++)
        if (std::abs(1-(x_m[i+1]-x_m[i])/(x_m[1]-x_m[0])) > 1e-9) {
            constantSpacing_m = false;
            return;
        }
    for (unsigned int i = 0; i < y_m.size()-1; i++)
        if (std::abs(1-(y_m[i+1]-y_m[i])/(y_m[1]-y_m[0])) > 1e-9) {
            constantSpacing_m = false;
            return;
        }
    for (unsigned int i = 0; i < z_m.size()-1; i++)
        if (std::abs(1-(z_m[i+1]-z_m[i])/(z_m[1]-z_m[0])) > 1e-9) {
            constantSpacing_m = false;
            return;
        }
}

Mesh::Iterator ThreeDGrid::getNearest(const double* position) const {
    std::vector<int> index(3);
    lowerBound(position[0], index[0],
               position[1], index[1],
               position[2], index[2]);
    if (index[0] < xSize_m-1 && index[0] >= 0)
        index[0] += (2*(position[0] - x_m[index[0]])  >
                     x_m[index[0]+1]-x_m[index[0]] ? 2 : 1);
    else
        index[0]++;
    if (index[1] < ySize_m-1 && index[1] >= 0)
        index[1] += (2*(position[1] - y_m[index[1]])  >
                    y_m[index[1]+1]-y_m[index[1]] ? 2 : 1);
    else
        index[1]++;
    if (index[2] < zSize_m-1 && index[2] >= 0)
        index[2] += (2*(position[2] - z_m[index[2]])  >
                    z_m[index[2]+1]-z_m[index[2]] ? 2 : 1);
    else
        index[2]++;
    if (index[0] < 1)
        index[0] = 1;
    if (index[1] < 1)
        index[1] = 1;
    if (index[2] < 1)
        index[2] = 1;
    if (index[0] > xSize_m)
        index[0] = xSize_m;
    if (index[1] > ySize_m)
        index[1] = ySize_m;
    if (index[2] > zSize_m)
        index[2] = zSize_m;
    return Mesh::Iterator(index, this);
}
}

