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

#ifndef _CLASSIC_FIELDS_NDGRID_HH_
#define _CLASSIC_FIELDS_NDGRID_HH_

#include <algorithm>
#include <cmath>
#include <ostream>
#include <vector>

#include "Fields/Interpolation/Mesh.h"

namespace interpolation {

/** NDGrid holds grid information for a rectilinear grid in some arbitrary
 *  dimensional space.
 *
 *  NDGrid holds the necessary information for regular and irregular grid data
 *  for interpolation routines (but it all has to be rectilinear).
 *
 *  Where an irregular grid is used, we use std::lower_bound to find mapping
 *  from a generic point to the grid points.
 *
 *  NDGrid holds a list of pointers to VectorMaps that use the grid for
 *  interpolation and functions to Add and Remove VectorMaps. If the Remove
 *  function removes the last VectorMap, then the ThreeDGrid is deleted. It's a
 *  sort of smart pointer.
 *
 *  In general NDGrid is not bound checked.
 */

class NDGrid : public Mesh {
  public:
    class Iterator;

    /** Inheritable copy constructor */
    inline NDGrid* clone();

    /** Not implemented */
    Mesh* dual() const;

    /** Build a default, empty grid with zero dimension */
    NDGrid();

    /** Build a regular grid
     *
     *  @param nDimensions: number of dimensions of the grid
     *  @param size: array of ints of length nDimensions. Gives the number of
     *               grid points in each dimension. Caller owns the memory
     *               pointed to by size. Must be at least 1 grid point in each
     *               dimension.
     *  @param spacing: array of doubles of length nDimensions. Gives the grid
     *               spacing in each dimension. Caller owns memory pointed to by
     *               spacing.
     *  @param min: array of doubles of length nDimensions. Gives the coordinate
     *              of the zeroth element (origin). Caller owns the memory
     *              pointed to by min.
     */
    NDGrid(int nDimensions, int* size, double* spacing, double* min);

    /** Build an irregular grid
     *
     *  @param size: The number of grid points in each dimension. This vector
     *               must be the same length as gridCoordinates vector.
     *  @param gridCoordinates: vector of doubles, with the same length as size.
     *               Set the coordinates of the grid in each dimension. Caller
     *               owns the memory pointed to by gridCoordinates.
     */
    NDGrid(std::vector<int> size, std::vector<const double *> gridCoordinates);

    /** Build an irregular grid
     *
     *  @param gridCoordinates: Vector of size equal to the dimension of the
     *                          grid. Element gridCoordinates[i][j] is the
     *                          position of the jth grid coordinate in the ith
     *                          dimension. i.e. to find the coordinates of the
     *                          grid point (23, 17, 21) we would lookup
     *                             gridCoordinate[0][23]
     *                             gridCoordinate[1][17]
     *                             gridCoordinate[2][21]
     */
    NDGrid(std::vector< std::vector<double> > gridCoordinates);

    /** Destructor */
    ~NDGrid() {;}

    /** Get the coordinate of a set of grid points
     *
     *  @param index: the index of the grid point to be found; indexed from 1.
     *         Should be in range 1 <= index <= length
     *  @param dimension: the dimension of the grid point to be found. Should be
     *         in range 0 <= dimension < grid dimension
     *
     *  Note the coord is NOT bound checked - so it is a memory error to go
     *  out of the bounds.
     */
    inline double& coord(const int& index, const int& dimension);

    /** Get the coordinate of a set of grid points
     *
     *  @param index: the index of the grid point to be found; indexed from 1.
     *         Should be in range 1 <= index <= length
     *  @param dimension: the dimension of the grid point to be found. Should be
     *         in range 0 <= dimension < grid dimension
     *
     *  Note the call is NOT bound checked - so it is a memory error to go
     *  out of the bounds.
     */
    inline const double& coord(const int& index, const int& dimension) const;

    /** Get the size of the grid along a given dimension
     *
     *  @param dimension: the dimension of the grid point to be found. Should be
     *         in range 0 <= dimension < grid dimension
     *
     *  Note the call is NOT bound checked - so it is a memory error to go
     *  out of the bounds.
     */
    inline int size(const int& dimension) const;

    /** Get the vector of grid points along a particular dimension
     *
     *  @param dimension: the dimension of the grid point to be found. Should be
     *         in range 0 <= dimension < grid dimension
     *
     *  Note the call is NOT bound checked - so it is a memory error to go
     *  out of the bounds.
     */
    inline std::vector<double> coordVector(const int& dimension) const;

    /** Get an array of grid points along a particular dimension
     *
     *  @param dimension: the dimension of the grid point to be found. Should be
     *         in range 0 <= dimension < grid dimension
     *
     *  Note the call is NOT bound checked - so it is a memory error to go
     *  out of the bounds.
     *  @returns an array of length this->size(dimension); caller owns the 
     *  returned memory.
     */
    double* newCoordArray  (const int& dimension) const;

    /** Get the grid point of the grid lower bound in a particular dimension
     *
     *  @param x: position in the grid coordinate
     *  @param dimension: the grid coordinate to search along
     *  @param xIndex: will be set with the resultant grid index
     *
     *  Sets xIndex to the index of the largest mesh position that is less than
     *  x in the given dimension. If x is less than all mesh positions, sets
     *  xIndex to -1.
     */
    inline void coordLowerBound(const double& x, const int& dimension, int& xIndex) const;

    /** Get the grid point of the grid lower bound
     *
     *  @param pos: position in the grid coordinate. Should have a length equal
     *         to the dimension of the grid.
     *  @param xIndex: will be set with the resultant grid index. Should have a
     *         lenght equal to the dimension of the grid.
     *
     *  Sets xIndex to the index of the largest mesh position that is less than
     *  pos in all dimension.
     *
     *  Note that dimension and xIndex are not bound checked.
     */
    inline void lowerBound(const std::vector<double>& pos, std::vector<int>& xIndex) const;

    /** Return the lowest grid coordinate in the given dimension */
    inline double  min(const int& dimension) const;

    /** Return the highest grid coordinate in the given dimension */
    inline double  max(const int& dimension) const;

    /** Reset the coordinates in the dimension to a new value
     *
     *  @param dimension: the dimension of the coordinates that will be reset
     *  @param nCoords: the length of array x
     *  @param x: array of points. Caller owns the memory allocated to x.
     */
    inline void setCoord(int dimension, int nCoords, double * x);

    /** Returns the origin of the mesh (lowest possible index in all dimensions)
     */
    inline Mesh::Iterator begin() const;

    /** Returns the end of the mesh (highest possible index in all dimensions)
     */
    inline Mesh::Iterator end() const;

    /** Get the position of the iterator
     *  
     *  @param it: iterator on this Mesh
     *  @param position: double array of length equal to the mesh dimension.
     *             Filled with the iterator position. Caller owns the memory
     *             to which position points.
     */
    inline void getPosition(const Mesh::Iterator& it, double * position) const;

    /** Get the dimension of the Mesh */
    inline int getPositionDimension() const;

    /** Returns true if the Mesh is a regular grid */
    inline bool getConstantSpacing() const;

    /** Set to true to flag the Mesh as a regular grid
     *
     *  @param spacing: flag indicating if the Mesh is regular.
     *
     *  Nb: If this is set for an irregular grid... I don't know what will
     *  happen
     */
    inline void setConstantSpacing(bool spacing);

    /** Automatically detector whether the Mesh is a regular grid
     *
     *  Loops over all grid points and checks the step size is regular, with
     *  tolerance given by tolernance_m
     */
    void setConstantSpacing(double tolerance_m = 1e-9);
 
    /** Return an integer corresponding to the iterator position
     *
     *  @param lhs: iterator whose position is returned.
     *  
     *  @returns an integer between 0 and number of grid points; 0 is
     *  this->begin(), and each call to iterator++ increments the integer by 1
     */
    int toInteger(const Mesh::Iterator& lhs) const;

    /** Return an iterator corresponding to nearest mesh point to position
     *
     *  @param position: Array of length dimension, giving a position near the
     *                   grid.
     *  
     *  @returns iterator that points to the grid location closest to position.
     */
    Mesh::Iterator getNearest(const double* position) const;

protected:

    //Change position
    virtual Mesh::Iterator& addEquals(Mesh::Iterator& lhs, int difference) const;
    virtual Mesh::Iterator& subEquals(Mesh::Iterator& lhs, int difference) const;
    virtual Mesh::Iterator& addEquals(Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;
    virtual Mesh::Iterator& subEquals(Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;
    virtual Mesh::Iterator& addOne   (Mesh::Iterator& lhs) const;
    virtual Mesh::Iterator& subOne   (Mesh::Iterator& lhs) const;
    //Check relative position
    virtual bool isGreater(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;

private:
    std::vector< std::vector<double> > coord_m;
    std::vector<VectorMap*>            maps_m;
    bool                               constantSpacing_m;

    friend Mesh::Iterator  operator++(Mesh::Iterator& lhs, int);
    friend Mesh::Iterator  operator--(Mesh::Iterator& lhs, int);
    friend Mesh::Iterator& operator++(Mesh::Iterator& lhs);
    friend Mesh::Iterator& operator--(Mesh::Iterator& lhs);
    friend Mesh::Iterator  operator- (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator  operator+ (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator& operator-=(Mesh::Iterator& lhs,       const Mesh::Iterator& rhs);
    friend Mesh::Iterator& operator+=(Mesh::Iterator& lhs,       const Mesh::Iterator& rhs);

    friend bool operator==(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator!=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator>=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator<=(const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator< (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend bool operator> (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
};

double& NDGrid::coord(const int& index, const int& dimension) {
    return coord_m[dimension][index-1];
}

const double& NDGrid::coord(const int& index, const int& dimension) const {
    return coord_m[dimension][index-1];
}

int NDGrid::size(const int& dimension) const {
    return coord_m[dimension].size();
}

std::vector<double> NDGrid::coordVector(const int& dimension) const {
    return coord_m[dimension];
}

void NDGrid::coordLowerBound(const double& x, const int& dimension, int& xIndex) const {
    if (constantSpacing_m) {
        double x0 = coord_m[dimension][0];
        double x1 = coord_m[dimension][1];
        xIndex = static_cast<int>(std::floor((x - x0)/(x1-x0)));
        int coordSize = static_cast<int>(coord_m[dimension].size()); 
        if (xIndex >= coordSize) {
            xIndex = coordSize-1;
        } else if (xIndex < -1) {
            xIndex = -1;
        }
    } else {
        const std::vector<double>& c_t(coord_m[dimension]);
        xIndex = std::lower_bound(c_t.begin(), c_t.end(), x)-c_t.begin()-1;
    }
}

void NDGrid::lowerBound(const std::vector<double>& pos, std::vector<int>& xIndex) const {
    xIndex = std::vector<int> (pos.size());
    for(unsigned int i=0; i < pos.size(); i++) {
        coordLowerBound(pos[i], i, xIndex[i]);
    }
}

double NDGrid::min(const int& dimension) const {
    return coord_m[dimension][0];
}

double NDGrid::max(const int& dimension) const {
    return coord_m[dimension][coord_m[dimension].size()-1];
}

NDGrid* NDGrid::clone() {
    return new NDGrid(*this);
}

void NDGrid::setCoord(int dimension, int nCoords, double * x) {
    coord_m[dimension] = std::vector<double>(x, x+nCoords);
}

Mesh::Iterator NDGrid::begin() const {
    return Mesh::Iterator(std::vector<int>(coord_m.size(), 1), this);
}

Mesh::Iterator NDGrid::end() const {
    if (coord_m.size() == 0) {
        return Mesh::Iterator(std::vector<int>(), this);
    }
    std::vector<int> end(coord_m.size(), 1);
    end[0] = coord_m[0].size()+1;
    return Mesh::Iterator(end, this);
}

void NDGrid::getPosition(const Mesh::Iterator& it, double * position) const {
    for (unsigned int i=0; i<it.getState().size(); i++)
        position[i] = coord_m[i][it[i]-1];
}

int NDGrid::getPositionDimension() const {
    return coord_m.size();
}

bool NDGrid::getConstantSpacing() const {
    return constantSpacing_m;
}

void NDGrid::setConstantSpacing(bool spacing) {
    constantSpacing_m = spacing;
}

} // namespace interpolation

#endif // _CLASSIC_FIELDS_NDGRID_HH_
