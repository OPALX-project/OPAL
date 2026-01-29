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

#ifndef _CLASSIC_FIELDS_INTERPOLATOR3DGRIDTO3D_HH_
#define _CLASSIC_FIELDS_INTERPOLATOR3DGRIDTO3D_HH_

#include "Fields/Interpolation/Mesh.h"
#include "Fields/Interpolation/VectorMap.h"
#include "Fields/Interpolation/Interpolator3dGridTo1d.h"
#include "Fields/Interpolation/TriLinearInterpolator.h"

namespace interpolation {

/** Interpolator3dGridTo3d interpolates from 3d grid to a 3d vector
 *
 *  Wraps three Interpolator3dGridTo1d, one for each variable in the output
 *  vector B. At the moment the wrapped Interpolator3dGridTo1d is encoded by
 *  means of an enumeration, with the only possible value trilinear interpolator
 *
 *  Function data for each of the 3d output data on the mesh is written as 3d
 *  arrays Bx, By, Bz and we make a TriLinear interpolator for each array
 *
 *  Could be that we make Interpolator3dGridTo3d an interface class and then
 *  make specific implementations for different interpolation routines (of which
 *  three trilinear interpolators is an implementation of the interface).
 */


class Interpolator3dGridTo3d : public VectorMap {
 public:
  /** Enumerator encoding possible interpolation routines */
  enum interpolationAlgorithm {TRILINEAR};

  /** Constructor for grids with constant spacing
   *
   *  \param grid 3d grid that indexes the data. *this is added to the smart
   *         pointer in the mesh
   *  \param Bx value[0] data. 3D arrays go like [index_x][index_y][index_z],
   *            and Interpolator3dGridTo3d now owns this memory
   *  \param By value[1] data. 3D arrays go like [index_x][index_y][index_z],
   *            and Interpolator3dGridTo3d now owns this memory
   *  \param Bz value[2] data. 3D arrays go like [index_x][index_y][index_z],
   *            and Interpolator3dGridTo3d now owns this memory
   *  \param algo interpolation algorithm (defaults to TRILINEAR)
   */
  inline Interpolator3dGridTo3d(ThreeDGrid* grid,
                         double ***Bx, double ***By, double ***Bz,
                         interpolationAlgorithm algo = TRILINEAR);

  /** Copy constructor deep copies 1d interpolators and mesh
   *
   *  Note this makes a whole bunch of extra meshes because each interpolator
   *  deep copy makes a new mesh - that's a bit of a mess
   */
  Interpolator3dGridTo3d(const Interpolator3dGridTo3d& interpolator);

  /** Delete member interpolators and remove *this from the mesh smart pointer
   */
  ~Interpolator3dGridTo3d() {clear();}

  /** Copy function (can be called on parent class) */
  Interpolator3dGridTo3d* clone() const
                                     {return new Interpolator3dGridTo3d(*this);}

  /** Return the interpolated data
   *
   *  Calls each child interpolator in turn and sets the value. Note that this
   *  means we get 3 sets of bound checks plus bound checking on the parent
   *  which is a bit unpleasant. I think that the actual interpolation however
   *  does have to be done independently for each variable.
   */
  void function(const double Point[3], double Value[3]) const;

  /** Call function at a particular point in the mesh */
  inline virtual void function
                            (const Mesh::Iterator& point, double* value) const {
     VectorMap::function(point, value);
  }

  /** Do not use (just raises exception) - der
   */
  void functionPrime(const double Point[3], double Value[3], int axis) const;

  /** Number of x coordinates in the grid */
  inline int getNumberOfXCoords() const;

  /** Number of y coordinates in the grid */
  inline int getNumberOfYCoords() const;

  /** Number of z coordinates in the grid */
  inline int getNumberOfZCoords() const;

  /** Dimension of input points */
  inline unsigned int getPointDimension()  const;

  /** Dimension of output values */
  inline unsigned int getValueDimension()  const;

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

  /** Set function and mesh data
   */
  void setAll(ThreeDGrid* grid, double *** Bx, double *** By, double *** Bz,
           interpolationAlgorithm algo = TRILINEAR);

  /** Clear all private data */
  inline void clear();

 protected:
  ThreeDGrid             *coordinates_m;
  Interpolator3dGridTo1d *interpolator_m[3];
};

int Interpolator3dGridTo3d::getNumberOfXCoords() const {
    return coordinates_m->xSize();
}

int Interpolator3dGridTo3d::getNumberOfYCoords() const {
    return coordinates_m->ySize();
}

int Interpolator3dGridTo3d::getNumberOfZCoords() const {
    return coordinates_m->zSize();
}

unsigned int Interpolator3dGridTo3d::getPointDimension()  const {
    return 3;
}

unsigned int Interpolator3dGridTo3d::getValueDimension()  const {
    return 3;
}

ThreeDGrid* Interpolator3dGridTo3d::getMesh() const {
    return coordinates_m;
}

void Interpolator3dGridTo3d::setGrid(ThreeDGrid* grid) {
    if (coordinates_m != NULL)
        coordinates_m->remove(this);
    grid->add(this);
    coordinates_m = grid;
}

void Interpolator3dGridTo3d::setX(int nCoords, double* x) {
    if (coordinates_m != NULL)
        coordinates_m->setX(nCoords, x);
}

void Interpolator3dGridTo3d::setY(int nCoords, double* y) {
    if (coordinates_m != NULL)
        coordinates_m->setY(nCoords, y);
}

void Interpolator3dGridTo3d::setZ(int nCoords, double* z) {
    if (coordinates_m != NULL)
        coordinates_m->setZ(nCoords, z);
}

inline Interpolator3dGridTo3d::Interpolator3dGridTo3d(ThreeDGrid* grid,
                       double ***Bx, double ***By, double ***Bz,
                       interpolationAlgorithm /*algo*/)
      : coordinates_m(NULL) {
    for (int i = 0; i < 3; i++)
        interpolator_m[i] = NULL;
    setAll(grid, Bx, By, Bz);
}

void Interpolator3dGridTo3d::clear() {
    for (int i = 0; i < 3; i++)
        delete interpolator_m[i];
    coordinates_m->remove(this);
}

}

#endif
