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

#ifndef _CLASSIC_FIELDS_THREEDGRID_HH_
#define _CLASSIC_FIELDS_THREEDGRID_HH_

#include <math.h>

#include <algorithm>
#include <vector>

#include "Fields/Interpolation/Mesh.h"

namespace interpolation {

/** ThreeDGrid holds grid information for a rectangular grid used in e.g.
 *  interpolation
 *
 *  ThreeDGrid holds the necessary information for regular and irregular grid
 *  data for interpolation routines.
 *
 *  Constructor is provided to generate a regular grid given grid spacing etc
 *  and an irregular grid given lists of x, y, z data for the grid points.
 *  Nearest either calculates nearest grid point a priori for regular grids or
 *  uses std::lower_bound to calculate nearest grid point (list bisection or
 *  somesuch). Controlled by _constantSpacing flag. The regular spacing lookup
 *  is quicker, so prefer constant spacing.
 *
 *  ThreeDGrid holds a list of pointers to VectorMaps that use the grid for
 *  interpolation and functions to Add and Remove VectorMaps. If the Remove
 *  function removes the last VectorMap, then the ThreeDGrid is deleted. It's a
 *  sort of smart pointer.
 */
class ThreeDGrid : public Mesh {
  public:
    class Iterator;

    /** Deep copy the grid */
    ThreeDGrid * clone() {return new ThreeDGrid(*this);}

    /** Not implemented (returns NULL) */
    Mesh * dual () const;

    /** Default constructor */
    ThreeDGrid();

    /** Make a regular rectangular ThreeDGrid
     *
     *  Set up a ThreeDGrid with regular grid spacing (same distance between
     *  every point)
     *  \param dX step size in x; should be positive
     *  \param dY step size in y; should be positive
     *  \param dZ step size in z; should be positive
     *  \param minX x position of grid lower edge
     *  \param minY y position of grid lower edge
     *  \param minZ z position of grid lower edge
     *  \param numberOfXCoords number of x coordinates on grid; should be > 1
     *  \param numberOfYCoords number of y coordinates on grid; should be > 1
     *  \param numberOfZCoords number of z coordinates on grid; should be > 1
     */
    ThreeDGrid(double dX, double dY, double dZ,
               double minX, double minY, double minZ,
               int numberOfXCoords, int numberOfYCoords, int numberOfZCoords);

    /** Make a irregular rectangular ThreeDGrid
     *
     *  Set up a ThreeDGrid with irregular grid spacing
     *  \param xSize number of points in x
     *  \param x x points array
     *  \param ySize number of points in y
     *  \param y y points array
     *  \param zSize number of points in Z
     *  \param z z point array
     *  ThreeDGrid deep copies the data pointed to by x, y, z (i.e. it does not
     *  take ownership of these arrays).
     */
    ThreeDGrid(int xSize, const double *x,
               int ySize, const double *y,
               int zSize, const double *z);

    /** Make a irregular rectangular ThreeDGrid
     *
     *  Set up a ThreeDGrid with irregular grid spacing
     *  \param x x points array
     *  \param y y points array
     *  \param z z point array
     */
    ThreeDGrid(std::vector<double> x,
               std::vector<double> y,
               std::vector<double> z);

    /** Destructor */
    ~ThreeDGrid();

    /** Get ith coordinate in x
     *
     *  Round bracket indexing starts at 1 goes to nCoords
     */
    inline double& x(const int& i) {return x_m[i-1];}

    /** Get ith coordinate in y
     *
     *  Round bracket indexing starts at 1 goes to nCoords
     */
    inline double& y(const int& j) {return y_m[j-1];}

    /** Get ith coordinate in z
     *
     *  Round bracket indexing starts at 1 goes to nCoords
     */
    inline double& z(const int& k) {return z_m[k-1];}

    /** Get ith coordinate in x
     *
     *  Round bracket indexing starts at 1 goes to nCoords
     */
    inline const double& x(const int& i) const {return x_m[i-1];}

    /** Get ith coordinate in y
     *
     *  Round bracket indexing starts at 1 goes to nCoords
     */
    inline const double& y(const int& j) const {return y_m[j-1];}

    /** Get ith coordinate in z
     *
     *  Round bracket indexing starts at 1 goes to nCoords
     */
    inline const double& z(const int& j) const {return z_m[j-1];}

    /** Get size of x data */
    inline int     xSize() const {return static_cast<int>(x_m.size());}

    /** Get size of y data */
    inline int     ySize() const {return static_cast<int>(y_m.size());}

    /** Get size of z data */
    inline int     zSize() const {return static_cast<int>(z_m.size());}

    /** Get vector containing x grid data */
    std::vector<double> xVector() {return std::vector<double>(x_m);}

    /** Get vector containing y grid data */
    std::vector<double> yVector() {return std::vector<double>(y_m);}

    /** Get vector containing z grid data */
    std::vector<double> zVector() {return std::vector<double>(z_m);}

    /** Allocate a new array containing x grid data */
    inline double* newXArray();

    /** Allocate a new array containing y grid data */
    inline double* newYArray();

    /** Allocate a new array containing z grid data */
    inline double* newZArray();

    /** Find the index of the nearest point less than x
     *
     *  Note that in the case of a regular grid there is no bound checking; in
     *  the case of an irregular grid we always return an element in the grid
     */
    inline void xLowerBound(const double& x, int& xIndex) const;

    /** Find the index of the nearest point less than y
     *
     *  Note that in the case of a regular grid there is no bound checking; in
     *  the case of an irregular grid we always return an element in the grid
     */
    inline void yLowerBound(const double& y, int& yIndex) const;

    /** Find the index of the nearest point less than z
     *
     *  Note that in the case of a regular grid there is no bound checking; in
     *  the case of an irregular grid we always return an element in the grid
     */
    inline void zLowerBound(const double& z, int& zIndex) const;

    /** Find the indices of the nearest point on the grid less than x, y, z
     *
     *  Note that in the case of a regular grid there is no bound checking; in
     *  the case of an irregular grid we always return an element in the grid
     */
    inline void lowerBound(const double& x, int& xIndex,
                           const double& y, int& yIndex,
                           const double& z, int& zIndex) const;

    /** Find the iterator of the nearest point on the grid less than x, y, z
     *
     *  Note that in the case of a regular grid there is no bound checking; in
     *  the case of an irregular grid we always return an element in the grid
     */
    inline void lowerBound(const double& x, const double& y, const double& z,
                           Mesh::Iterator& it) const;

    /** Return the smallest value of x in the grid */
    inline double minX() const;

    /** Return the largest value of x in the grid */
    inline double maxX() const;

    /** Return the smallest value of y in the grid */
    inline double minY() const;

    /** Return the largest value of y in the grid */
    inline double maxY() const;

    /** Return the smallest value of z in the grid */
    inline double minZ() const;

    /** Return the largest value of z in the grid */
    inline double maxZ() const;

    /** Reset x grid points - note may need to set SetConstantSpacing() */
    inline void setX(int nXCoords, double * x);

    /** Reset y grid points - note may need to set SetConstantSpacing() */
    inline void setY(int nYCoords, double * y);

    /** Reset z grid points - note may need to set SetConstantSpacing() */
    inline void setZ(int nZCoords, double * z);

    /** Add *map to the maps_m list if it has not already been added */
    void add(VectorMap* map);

    /** Remove *map from the maps_m list if it has not already been removed. If
     *  there are no more maps in the list, delete this*/
    void remove(VectorMap* map);

    /** Return an iterator pointing at the first element in the ThreeDGrid */
    Mesh::Iterator begin() const;

    /** Return an iterator pointing at the last+1 element in the ThreeDGrid */
    Mesh::Iterator end()   const;

    /** Fill position with the position of it
     *
     *  \param it iterator pointing at a grid position
     *  \param position array of length >= 3 that will be filled with the
     *         position of it
     */
    virtual void getPosition(const Mesh::Iterator& it, double * position) const;

    /** Returns 3, the dimension of the space of the grid */
    inline int getPositionDimension() const;

    /** Converts from iterator to integer value
     *
     *  0 corresponds to the mesh beginning. integers index (right most index)
     *  indexes z, middle index is y, most significant (left most index) indexes
     *  x.
     */
    inline int toInteger(const Mesh::Iterator& lhs) const;

    /** Set to true to use constant spacing
     *
     *  If constant spacing is true, assume regular grid with spacing given by
     *  separation of first two elements in each dimension
     */
    inline void setConstantSpacing(bool spacing);

    /** Autodetect constant spacing with float tolerance of 1e-9
     */
    void setConstantSpacing();

    /** Return true if using constant spacing */
    inline bool getConstantSpacing() const;

    /** Return the point in the mesh nearest to position */
    Mesh::Iterator getNearest(const double* position) const;

    /** Custom LowerBound routine like std::lower_bound
     *
     *  Make a binary search to find the element in vec with vec[i] <= x.
     *  \param vec STL vector object (could make this any sorted iterator...)
     *  \param x object for which to search.
     *  \param index vectorLowerBound sets index to the position of the element. If
     *    x < vec[0], vectorLowerBound fills with -1.
     */
    static void vectorLowerBound(std::vector<double> vec, double x, int& index);

  protected:
    // Change position
    virtual Mesh::Iterator& addEquals(Mesh::Iterator& lhs,
                                      int difference) const;
    virtual Mesh::Iterator& subEquals(Mesh::Iterator& lhs,
                                      int difference) const;
    virtual Mesh::Iterator& addEquals
                        (Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;
    virtual Mesh::Iterator& subEquals
                        (Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;
    virtual Mesh::Iterator& addOne(Mesh::Iterator& lhs) const;
    virtual Mesh::Iterator& subOne(Mesh::Iterator& lhs) const;
    // Check relative position
    virtual bool isGreater
                  (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs) const;

  private:
    std::vector<double>     x_m;
    std::vector<double>     y_m;
    std::vector<double>     z_m;
    int                     xSize_m;
    int                     ySize_m;
    int                     zSize_m;
    std::vector<VectorMap*> maps_m;
    bool                    constantSpacing_m;

    friend Mesh::Iterator  operator++(Mesh::Iterator& lhs, int);
    friend Mesh::Iterator  operator--(Mesh::Iterator& lhs, int);
    friend Mesh::Iterator& operator++(Mesh::Iterator& lhs);
    friend Mesh::Iterator& operator--(Mesh::Iterator& lhs);
    friend Mesh::Iterator  operator-
                        (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator  operator+
                        (const Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator& operator-=
                        (Mesh::Iterator& lhs, const Mesh::Iterator& rhs);
    friend Mesh::Iterator& operator+=
                        (Mesh::Iterator& lhs, const Mesh::Iterator& rhs);

    friend bool operator==(const Mesh::Iterator& lhs,
                           const Mesh::Iterator& rhs);
    friend bool operator!=(const Mesh::Iterator& lhs,
                           const Mesh::Iterator& rhs);
    friend bool operator>=(const Mesh::Iterator& lhs,
                           const Mesh::Iterator& rhs);
    friend bool operator<=(const Mesh::Iterator& lhs,
                           const Mesh::Iterator& rhs);
    friend bool operator< (const Mesh::Iterator& lhs,
                           const Mesh::Iterator& rhs);
    friend bool operator> (const Mesh::Iterator& lhs,
                           const Mesh::Iterator& rhs);
};

double* ThreeDGrid::newXArray() {
    double *x = new double[x_m.size()];
    for (unsigned int i = 0; i < x_m.size(); i++)
        x[i] = x_m[i];
    return x;
}

double* ThreeDGrid::newYArray() {
    double *y = new double[y_m.size()];
    for (unsigned int i = 0; i < y_m.size(); i++)
        y[i] = y_m[i];
    return y;
}

double* ThreeDGrid::newZArray() {
    double *z = new double[z_m.size()];
    for (unsigned int i = 0; i < z_m.size(); i++)
        z[i] = z_m[i];
    return z;
}

void ThreeDGrid::lowerBound(const double& x, int& xIndex,
                            const double& y, int& yIndex,
                            const double& z, int& zIndex) const {
    xLowerBound(x, xIndex);
    yLowerBound(y, yIndex);
    zLowerBound(z, zIndex);
}

void ThreeDGrid::lowerBound(const double& x,
                            const double& y,
                            const double& z,
                            Mesh::Iterator& it) const {
    xLowerBound(x, it[0]);
    yLowerBound(y, it[1]);
    zLowerBound(z, it[2]);
    it[0]++;
    it[1]++;
    it[2]++;
}

void ThreeDGrid::xLowerBound(const double& x, int& xIndex) const {
    if (constantSpacing_m)
        xIndex = static_cast<int>(std::floor((x - x_m[0])/(x_m[1]-x_m[0]) ));
    else {
        vectorLowerBound(x_m, x, xIndex);
    }
}


void ThreeDGrid::yLowerBound(const double& y, int& yIndex) const {
    if (constantSpacing_m)
        yIndex = static_cast<int>(std::floor((y - y_m[0])/(y_m[1]-y_m[0])));
    else
        vectorLowerBound(y_m, y, yIndex);
}

void ThreeDGrid::zLowerBound(const double& z, int& zIndex) const {
    if (constantSpacing_m)
        zIndex = static_cast<int>(std::floor((z - z_m[0])/(z_m[1]-z_m[0])));
    else
        vectorLowerBound(z_m, z, zIndex);
}

double ThreeDGrid::minX() const {
    return x_m[0];
}

double ThreeDGrid::maxX() const {
    return x_m[xSize_m-1];
}

double ThreeDGrid::minY() const {
    return y_m[0];
}

double ThreeDGrid::maxY() const {
    return y_m[ySize_m-1];
}

double ThreeDGrid::minZ() const {
    return z_m[0];
}

double ThreeDGrid::maxZ() const {
    return z_m[zSize_m-1];
}

void ThreeDGrid::setX(int nXCoords, double * x) {
    x_m = std::vector<double>(x, x+nXCoords);
}

void ThreeDGrid::setY(int nYCoords, double * y) {
    y_m = std::vector<double>(y, y+nYCoords);
}

void ThreeDGrid::setZ(int nZCoords, double * z) {
    z_m = std::vector<double>(z, z+nZCoords);
}

int ThreeDGrid::getPositionDimension() const {
    return 3;
}

int ThreeDGrid::toInteger(const Mesh::Iterator& lhs) const {
    return (lhs.getState()[0]-1)*zSize_m*ySize_m+
           (lhs.getState()[1]-1)*zSize_m+
           (lhs.getState()[2]-1);
}

void ThreeDGrid::setConstantSpacing(bool spacing) {
    constantSpacing_m = spacing;
}

bool ThreeDGrid::getConstantSpacing() const {
    return constantSpacing_m;
}
}
#endif

