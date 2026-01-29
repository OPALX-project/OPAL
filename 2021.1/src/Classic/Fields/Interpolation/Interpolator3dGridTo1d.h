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

#ifndef _CLASSIC_FIELDS_INTERPOLATOR3DGRIDTO1D_H_
#define _CLASSIC_FIELDS_INTERPOLATOR3DGRIDTO1D_H_

#include "Fields/Interpolation/VectorMap.h"
#include "Fields/Interpolation/ThreeDGrid.h"

namespace interpolation {

/** Interpolator3dGridTo1d is an abstraction for lookup on a 3D mesh to get a 1D
 *  value.
 *
 *  Interpolator3dGridTo1d abstracts the lookup from a 3D mesh to get a 1D value
 *  This is a building block for higher dimension (e.g. field vector lookup)
 *  where we can do multiple lookups one for each field, at the cost of doing
 *  multiple access/interpolation operations to the grid.
 */
class Interpolator3dGridTo1d : public VectorMap {
 public:
  /** Constructor for grids with constant spacing
   *
   *  \param grid 3d mesh on which data is stored - adds a reference to *this
   *              into the mesh smart pointer thing.
   *  \param F function data with points on each element of the grid. Indexing
   *           goes like [index_x][index_y][index_z]. Interpolator3dGridTo1d now
   *           owns this memory
   */
  inline Interpolator3dGridTo1d(ThreeDGrid* grid, double ***F);

  /** Default constructor - initialises data to NULL
   *
   *  Note I have no NULL checking on most member functions so be careful
   */
  inline Interpolator3dGridTo1d();

  /** Deletes F data and removes *this from ThreeDGrid smart pointer thing. */
  inline virtual ~Interpolator3dGridTo1d();

  /** Call function at a particular point in the mesh */
  inline virtual void function
                            (const Mesh::Iterator& point, double* value) const {
     VectorMap::function(point, value);
  }

  /** Get the interpolated function data */
  virtual void function(const double Point[3], double Value[1]) const = 0;

  /** Copy constructor */
  virtual Interpolator3dGridTo1d* clone() const = 0;

  /** Number of x coordinates in the grid */
  inline int getNumberOfXCoords() const;

  /** Number of y coordinates in the grid */
  inline int getNumberOfYCoords() const;

  /** Number of z coordinates in the grid
   */
  inline int getNumberOfZCoords() const;

  /** Dimension of input points */
  inline unsigned int getPointDimension() const;

  /** Dimension of output values */
  inline unsigned int getValueDimension() const;

  /** Return a pointer to the mesh*/
  inline ThreeDGrid* getMesh() const;

  /** Reset the mesh
   *
   *  Removes reference to this from mesh smart pointer if appropriate
   */
  inline void setGrid(ThreeDGrid* grid);

  /** Set x coordinates in the mesh to an arbitrary set of points
   *
   *  If mesh is NULL, does nothing
   */
  inline void setX(int nCoords, double* x);

  /** Set y coordinates in the mesh to an arbitrary set of points
   *
   *  If mesh is NULL, does nothing
   */
  inline void setY(int nCoords, double* y);

  /** Set z coordinates in the mesh to an arbitrary set of points
   *
   *  If mesh is NULL, does nothing
   */
  inline void setZ(int nCoords, double* z);

  /** Set function data, deleting any existing function data
   */
  inline void setF(double*** inF);

  /** Delete any existing function data
   */
  void deleteFunc(double*** func);

  /** Return a pointer to the function data
   */
  inline double*** function() const;

  /** Set function and mesh data
   */
  inline void setAll(ThreeDGrid* grid, double *** F);

  /** Clear all private data */
  inline void clear();
 protected:
  ThreeDGrid   *coordinates_m;
  double     ***f_m;
};

Interpolator3dGridTo1d::Interpolator3dGridTo1d(ThreeDGrid* grid, double ***F)
      : coordinates_m(NULL), f_m(NULL) {
    setAll(grid, F);
}

Interpolator3dGridTo1d::Interpolator3dGridTo1d()
      : coordinates_m(NULL), f_m(NULL) {
}

Interpolator3dGridTo1d::~Interpolator3dGridTo1d() {
    clear();
}

int Interpolator3dGridTo1d::getNumberOfXCoords() const {
    return coordinates_m->xSize();
}

int Interpolator3dGridTo1d::getNumberOfYCoords() const {
    return coordinates_m->ySize();
}

int Interpolator3dGridTo1d::getNumberOfZCoords() const {
    return coordinates_m->zSize();
}

unsigned int Interpolator3dGridTo1d::getPointDimension()  const {
    return 3;
}

unsigned int Interpolator3dGridTo1d::getValueDimension()  const {
    return 1;
}

ThreeDGrid* Interpolator3dGridTo1d::getMesh() const {
    return coordinates_m;
}

void Interpolator3dGridTo1d::setF(double*** inF) {
    deleteFunc(f_m);
    f_m = inF;
}

void Interpolator3dGridTo1d::setGrid(ThreeDGrid* grid) {
    if (coordinates_m != NULL)
        coordinates_m->remove(this);
    grid->add(this);
    coordinates_m = grid;
}

void Interpolator3dGridTo1d::setX(int nCoords, double* x) {
    if (coordinates_m != NULL)
        coordinates_m->setX(nCoords, x);
}

void Interpolator3dGridTo1d::setY(int nCoords, double* y) {
    if (coordinates_m != NULL)
        coordinates_m->setY(nCoords, y);
}

void Interpolator3dGridTo1d::setZ(int nCoords, double* z) {
    if (coordinates_m != NULL)
        coordinates_m->setZ(nCoords, z);
}

double*** Interpolator3dGridTo1d::function() const {
    return f_m;
}

void Interpolator3dGridTo1d::setAll(ThreeDGrid* grid, double *** F) {
    setGrid(grid);
    setF(F);
}

void Interpolator3dGridTo1d::clear() {
    deleteFunc(f_m);
    coordinates_m->remove(this);
}

}
#endif  // _CLASSIC_FIELDS_INTERPOLATOR3DGRIDTO1D_HH_
