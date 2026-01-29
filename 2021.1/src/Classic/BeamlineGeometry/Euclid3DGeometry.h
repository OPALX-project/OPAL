/* 
 *  Copyright (c) 2014, Chris Rogers
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

#ifndef CLASSIC_BEAMLINEGEOMETRY_Euclid3DGeometry_HH
#define CLASSIC_BEAMLINEGEOMETRY_Euclid3DGeometry_HH

#include "BeamlineGeometry/Geometry.h"
#include "BeamlineGeometry/Euclid3D.h"

// Euclid3D represents an aribitrary 3-d rotation and displacement.
class Euclid3D;

/** class Euclid3DGeometry
 *
 *  Euclid3DGeometry is a straight line from input to output followed by a 
 *  rotation.
 */

class Euclid3DGeometry : public BGeometryBase {
  public:

    Euclid3DGeometry(Euclid3D transformation);
    Euclid3DGeometry(const Euclid3DGeometry &right);
    virtual ~Euclid3DGeometry();
    const Euclid3DGeometry &operator=(const Euclid3DGeometry &right);

    /// Get arc length.
    //  Return the length of the geometry, measured along the design arc.
    virtual double getArcLength() const;

    /// Get geometry length.
    //  Return or the design length of the geometry.
    //  Depending on the element this may be the arc length or the
    //  straight length.
    virtual double getElementLength() const;

    /// Set geometry length.
    //  Assign the design length of the geometry.
    //  Depending on the element this may be the arc length or the
    //  straight length.
    virtual void setElementLength(double length);

    /// Not Implemented - raises GeneralClassicException
    /// Get transform.
    //  Return the transform of the local coordinate system from the
    //  position [b]fromS[/b] to the position [b]toS[/b].
    virtual Euclid3D getTransform(double fromS, double toS) const;

    /// Get total transform from beginning to end
    //  Corresponds to the Euclid3D
    virtual Euclid3D getTotalTransform() const;
  private:
    Euclid3D transformation_m;
};

#endif // #ifndef CLASSIC_BEAMLINEGEOMETRY_Euclid3DGeometry_HH
