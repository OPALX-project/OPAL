/*
 *  Copyright (c) 2018, Martin Duy Tat
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

#include <cmath>
#include <vector>
#include "BeamlineGeometry/VarRadiusGeometry.h"
#include "BeamlineGeometry/Euclid3D.h"
#include "AbsBeamline/MultipoleTFunctions/CoordinateTransform.h"

Euclid3D VarRadiusGeometry::getTransform(double fromS, double toS) const {
    Euclid3D v;
    coordinatetransform::CoordinateTransform t(0.0, 0.0, 0.0,
                                               s_0_m, lambda_left_m,
                                               lambda_right_m, rho_m);
    double phifrom = acos(t.getUnitTangentVector(fromS)[1]);
    double phito = acos(t.getUnitTangentVector(toS)[1]);
    std::vector<double> ref_from = t.calcReferenceTrajectory(fromS);
    std::vector<double> ref_to = t.calcReferenceTrajectory(toS);
    v = Euclid3D::YRotation(-(phifrom + phito));
    v.setX(ref_to[0] - ref_from[0]);
    v.setZ(ref_to[1] - ref_from[1]);
    return v;
}
