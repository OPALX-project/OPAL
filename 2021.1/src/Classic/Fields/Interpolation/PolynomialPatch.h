/* 
 *  Copyright (c) 2015, Chris Rogers
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

#include "Fields/Interpolation/SquarePolynomialVector.h"
#include "Fields/Interpolation/VectorMap.h"
#include "Fields/Interpolation/Mesh.h"

#ifndef PolynomialPatch_hh
#define PolynomialPatch_hh

namespace interpolation {

/** \class PolynomialPatch
 *  \brief Patches together many SquarePolynomialVectors to
 *  make a multidimensional polynomial spline.
 *
 *  SquarePolynomialVectors are distributed on a mesh. 
 *  \param grid_points_ The grid points on which patches are defined.
 *         SquarePolynomialVector owns this memory (deletes on destructor).
 *  \param validity_region_ Defines whether a region is valid or not. This is
 *         defined as a mesh for compatibility with the
 *         SectorMagneticFieldMap (but possibly this is the wrong way to do
 *         it - better as a lookup, like bool is_valid(point), which could be
 *         inherited from the VectorMap class)
 *         SquarePolynomialVector owns this memory (deletes on destructor).
 *  \param points_ The PolynomialVector data
 *         SquarePolynomialVector owns this memory (deletes on destructor).
 *  \param point_dimension_ Dimension of the ordinate
 *  \param value_dimension_ Dimension of the abscissa
 */
class PolynomialPatch : public VectorMap {
  public:
    /** Construct the PolynomialPatch
     *
     *  \param grid_points_ the mesh over which the polynomials are
     *         distributed. PolynomialPatch takes ownership of this memory.
     *  \param validity_region the mesh over which the grid points are valid;
     *         the validity region is defined by min/max values of the grid.
     *         PolynomialPatch takes ownership of this memory.
     *  \param polynomials_ the set of polynomials which make up the patch.
     *         PolynomialPatch takes ownership of this memory.
     */
    PolynomialPatch(Mesh* grid_points_,
                    Mesh* validity_region,
                    std::vector<SquarePolynomialVector*> polynomials_);

    /** Default constructor leaves validity_region_ and grid_points_ as NULL
     *
     *  Note it is not safe to call some member functions e.g. function(...) on
     *  default constructed PolynomialPatch.
     */
    PolynomialPatch();

    /** Destructor delets validity_region_, grid_points_ and points_ */
    ~PolynomialPatch();

    /** Deep copy the PolynomialPatch
     *
     *  This is the polymorphic copy constructor. Always returns a
     *  PolynomialPatch.
     */
    VectorMap* clone() const; //copy function

    /** Get the value at a given point 
     *
     *  \param point array of length point_dimension_. Caller owns this memory.
     *        Not bound checked - it wants to be fast.
     *  \param value array of length value_dimension_. Data gets overwritten by
     *        PolynomialPatch. Caller owns this memory. Not bound checked.
     */
    virtual void function(const double* point, double* value) const;

    /** Get the point dimension (length of the ordinate) */
    inline unsigned int getPointDimension() const {return point_dimension_;}

    /** Get the value dimension (length of the abscissa) */
    inline unsigned int getValueDimension() const {return value_dimension_;}

    /** Get the nearest SquarePolynomialVector to point
     *
     *  \param point array of length point_dimension_. Caller owns this memory.
     *        Not bound checked.
     */
    SquarePolynomialVector* getPolynomialVector(const double* point) const;

    /** Get the validity region mesh
     *
     *  Note that the PolynomialVectors are distributed on the dual of this mesh
     */
    Mesh* getMesh() const {return validity_region_;}

    /** Get all polynomials in the PolynomialPatch
     *
     *  This is a borrowed reference - PolynomialPatch still owns this memory.
     */
    std::vector<SquarePolynomialVector*> getPolynomials() const {return points_;}

  private:
    Mesh* validity_region_;
    Mesh* grid_points_;
    std::vector<SquarePolynomialVector*> points_;
    unsigned int point_dimension_;
    unsigned int value_dimension_;

    PolynomialPatch(const PolynomialPatch& rhs);
    PolynomialPatch& operator=(const PolynomialPatch& rhs);
};

}

#endif // PolynomialPatch_hh


