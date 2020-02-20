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

#ifndef _CLASSIC_FIELDS_VECTORMAP_HH_
#define _CLASSIC_FIELDS_VECTORMAP_HH_

#include <vector>

#include "Fields/Interpolation/ThreeDGrid.h"

namespace interpolation {

/** VectorMap is an abstract class that defines mapping from one vector to
 *  another.
 *
 *  VectorMap is base class primarily envisaged for mesh interpolation that
 *  defines interfaces to get interpolated data at a given point, access the
 *  mesh, etc.
 *
 *  Input and output dimensions (vector lengths) are arbitrary.
 */
class VectorMap {
 public:
    /** Pure virtual function to fill the array value with data evaluated at
     *  point.
     */
    virtual void function(const double* point, double* value) const = 0;

    /** Fill the array value with function data at a point on the mesh.
     *
     *  Implemented for default case where Mesh point dimension is the same as
     *  this point dimension (usual case). No default checking for array sizes.
     */
    inline virtual void function
                             (const Mesh::Iterator& point, double* value) const;

    /** Calculate F, appending output values to value_vec.
     *
     *  For each item in point_vec not in value_vec, calculate value_vec (urgh)
     */
    inline virtual void functionAppend(
                           const std::vector< std::vector<double> >& point_vec,
                           std::vector< std::vector<double> >& value_vec
                      ) const;

    /** Return true if point.size() is the same as this->PointDimension() */
    inline virtual bool  checkPoint(const std::vector<double>& point) const;

    /** Return true if value.size() is the same as this->ValueDimension() */
    inline virtual bool  checkValue(const std::vector<double>& value) const;

    /** Return the dimension of the point (input) */
    virtual unsigned int getPointDimension() const = 0;
    // would like to make static - but can't inherit static functions

    /** Return the dimension of the value (output) */
    virtual unsigned int getValueDimension() const = 0;
    // would like to make static - but can't inherit static functions

    /** Clone() is like a copy constructor - but copies the child class */
    virtual              VectorMap* clone() const = 0;

    /** Destructor */
    virtual             ~VectorMap() {;}

    /** Return the mesh used by the vector map or NULL if no mesh */
    virtual Mesh* getMesh() const {return NULL;}
  private:
};

bool VectorMap::checkPoint(const std::vector<double>& point) const {
  return (point.size() == this->getPointDimension());
}

bool VectorMap::checkValue(const std::vector<double>& value) const {
  return (value.size() == this->getValueDimension());
}

void VectorMap::function(const Mesh::Iterator& point, double* value) const {
    double PointA[this->getPointDimension()];
    point.getPosition(PointA);
    function(PointA, value);
}

void VectorMap::functionAppend
                       (const std::vector< std::vector<double> >& point_vec,
                        std::vector< std::vector<double> >& value_vec) const {
  for (size_t i = value_vec.size(); i < point_vec.size(); i++) {
    value_vec.push_back(std::vector<double>(getValueDimension()));
    function(&point_vec[i][0], &value_vec[i][0]);
  }
}
}
#endif  // _CLASSIC_FIELDS_VECTORMAP_HH_
// vi: set et ts=4 sw=4 sts=4:
// Local Variables:
// mode:c
// c-basic-offset: 4
// indent-tabs-mode: nil
// require-final-newline: nil
// End:
