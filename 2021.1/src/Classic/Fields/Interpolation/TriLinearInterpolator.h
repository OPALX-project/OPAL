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

#ifndef _SRC_LEGACY_INTERFACE_INTERPOLATION_TRILINEARINTERPOLATOR_HH_
#define _SRC_LEGACY_INTERFACE_INTERPOLATION_TRILINEARINTERPOLATOR_HH_

#include "Fields/Interpolation/Interpolator3dGridTo1d.h"
namespace interpolation {

/** TriLinearInterpolator performs a linear interpolation in x then y then z
 *
 *  Performs a linear interpolation in x then y then z returning a 1D value.
 */
class TriLinearInterpolator : public Interpolator3dGridTo1d {
  public:
    /** Constructor for grids with constant spacing
     *
     *  \param grid 3d mesh on which data is stored - adds a reference to *this
     *              into the mesh smart pointer thing.
     *  \param F function data with points on each element of the grid. Indexing
     *           goes like [index_x][index_y][index_z]. Interpolator3dGridTo1d now
     *           owns this memory
     *  All the data handling is done at Interpolator3dGridTo1d
     */
    inline TriLinearInterpolator(ThreeDGrid *grid, double ***F);

    /** Copy constructor
     *
     *  Deep copies the mesh and the function data
     */
    TriLinearInterpolator(const TriLinearInterpolator& tli);

    /** Destructor - removes reference from the mesh and from the function data
     */
    inline ~TriLinearInterpolator();

    /** Get the interpolated value of the function at some point
     *
     *  First does bound checking, then makes linear interpolations using the
     *  standard 1d interpolation formula \n
     *  \f$y(x) \approx \frac{\Delta y}{\Delta x} dx + y_0\f$ \n
     *  \f$y(x) \approx \frac{y_1(x_1)-y_0(x_0)}{dx}(x-x_0) + y_0\f$ \n
     *  Interpolate along 4 x grid lines to make a 2D problem, then interpolate
     *  along 2 y grid lines to make a 1D problem, then finally interpolate in z
     *  to get the value.
     */
    void function(const double Point[3], double Value[1]) const;

    /** Call function at a particular point in the mesh */
    inline virtual void function
                            (const Mesh::Iterator& point, double* value) const {
        VectorMap::function(point, value);
    }

    /** Copy function (can be called on parent class) */
    inline TriLinearInterpolator* clone() const;
};

TriLinearInterpolator::TriLinearInterpolator(ThreeDGrid *grid, double ***F)
    : Interpolator3dGridTo1d(grid, F) {
}

TriLinearInterpolator::~TriLinearInterpolator() {
}

TriLinearInterpolator* TriLinearInterpolator::clone() const {
    return new TriLinearInterpolator(*this);
}
}
#endif


