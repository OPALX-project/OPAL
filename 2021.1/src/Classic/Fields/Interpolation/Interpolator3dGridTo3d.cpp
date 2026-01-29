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

#include "Fields/Interpolation/Interpolator3dGridTo3d.h"
#include "Utilities/LogicalError.h"

namespace interpolation {

Interpolator3dGridTo3d::Interpolator3dGridTo3d
                                           (const Interpolator3dGridTo3d& rhs) {
    coordinates_m = new ThreeDGrid(*rhs.coordinates_m);
    for (int i = 0; i < 3; i++)
        interpolator_m[i] = rhs.interpolator_m[i]->clone();
}

void Interpolator3dGridTo3d::function
                                (const double Point[3], double Value[3]) const {
    if (Point[0] > coordinates_m->maxX() || Point[0] < coordinates_m->minX() ||
        Point[1] > coordinates_m->maxY() || Point[1] < coordinates_m->minY() ||
        Point[2] > coordinates_m->maxZ() || Point[2] < coordinates_m->minZ() ) {
        Value[0] = 0;
        Value[1] = 0;
        Value[2] = 0;
        return;
    }

    interpolator_m[0]->function(Point, &Value[0]);
    interpolator_m[1]->function(Point, &Value[1]);
    interpolator_m[2]->function(Point, &Value[2]);
}

void Interpolator3dGridTo3d::setAll(ThreeDGrid* grid,
                                 double *** Bx, double *** By, double *** Bz,
                                 interpolationAlgorithm algo) {
    if (coordinates_m != NULL)
      coordinates_m->remove(this);
    grid->add(this);
    coordinates_m = grid;

    for (int i = 0; i < 3; i++)
        if (interpolator_m[i] != NULL)
        delete interpolator_m[i];

    switch (algo) {
        case TRILINEAR:
            interpolator_m[0] = new TriLinearInterpolator(grid, Bx);
            interpolator_m[1] = new TriLinearInterpolator(grid, By);
            interpolator_m[2] = new TriLinearInterpolator(grid, Bz);
            break;
        default:
            throw(LogicalError(
                    "Interpolator3dGridTo3d::setAll",
                    "Did not recognise interpolation algorithm"
                   )
      );
    }
}
}
