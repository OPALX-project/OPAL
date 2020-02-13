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

#include "BeamlineGeometry/Euclid3DGeometry.h"

#include "Utilities/GeneralClassicException.h"
#include <cmath>

Euclid3DGeometry::Euclid3DGeometry(Euclid3D transformation)
    : BGeometryBase(), transformation_m(transformation) {
}

Euclid3DGeometry::Euclid3DGeometry(const Euclid3DGeometry &rhs)
    : BGeometryBase(), transformation_m(rhs.transformation_m) {
}

Euclid3DGeometry::~Euclid3DGeometry() {}

const Euclid3DGeometry &Euclid3DGeometry::operator=(const Euclid3DGeometry &rhs) {
    if (&rhs != this) {
        transformation_m = rhs.transformation_m;
    }
    return *this;
}

double Euclid3DGeometry::getArcLength() const {
    return sqrt(dot(transformation_m.getVector(),
                    transformation_m.getVector()));
}

double Euclid3DGeometry::getElementLength() const {
    return getArcLength();
}

void Euclid3DGeometry::setElementLength(double length) {
    Vector3D newVector = transformation_m.getVector()*(length/getArcLength());
    transformation_m.setDisplacement(newVector);
}

Euclid3D Euclid3DGeometry::getTransform(double /*fromS*/, double /*toS*/) const {
    throw GeneralClassicException("Euclid3DGeometry::getTransform", "Not implemented");
}

Euclid3D Euclid3DGeometry::getTotalTransform() const {
    return transformation_m;
}
