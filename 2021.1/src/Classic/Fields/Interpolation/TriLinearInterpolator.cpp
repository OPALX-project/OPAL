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

#include "Fields/Interpolation/TriLinearInterpolator.h"

namespace interpolation {

TriLinearInterpolator::TriLinearInterpolator(const TriLinearInterpolator& lhs):
  Interpolator3dGridTo1d(lhs)
{
    coordinates_m = new ThreeDGrid(*lhs.coordinates_m);
    f_m = new double**[coordinates_m->xSize()];
    for (int i = 0; i < coordinates_m->xSize(); i++) {
        f_m[i] = new double*[coordinates_m->ySize()];
        for (int j = 0; j < coordinates_m->ySize(); j++) {
            f_m[i][j] = new double[coordinates_m->zSize()];
            for (int k = 0; k < coordinates_m->zSize(); k++)
                f_m[i][j][k] = lhs.f_m[i][j][k];
        }
    }
    coordinates_m->add(this);
}

void TriLinearInterpolator::function
                                (const double Point[3], double Value[1]) const {
    int i, j, k;  // position indices
    coordinates_m->lowerBound(Point[0], i, Point[1], j, Point[2], k);
    // bound checking
    if (i + 2 > coordinates_m->xSize() ||
        j + 2 > coordinates_m->ySize() ||
        k + 2 > coordinates_m->zSize()) {
        Value[0] = 0;
        return;
    }
    if (i  < 0 || j < 0 || k < 0) {
        Value[0] = 0;
        return;
    }
    // interpolation in x
    double dx = (Point[0]-coordinates_m->x(i+1))/
                (coordinates_m->x(i+2)-coordinates_m->x(i+1));
    double f_x[2][2];
    f_x[0][0] = (f_m[i+1][j]  [k]   - f_m[i][j]  [k])   * dx + f_m[i][j]  [k];
    f_x[1][0] = (f_m[i+1][j+1][k]   - f_m[i][j+1][k])   * dx + f_m[i][j+1][k];
    f_x[0][1] = (f_m[i+1][j]  [k+1] - f_m[i][j]  [k+1]) * dx + f_m[i][j]  [k+1];
    f_x[1][1] = (f_m[i+1][j+1][k+1] - f_m[i][j+1][k+1]) * dx + f_m[i][j+1][k+1];
    // interpolation in y
    double f_xy[2];
    double dy = (Point[1]-coordinates_m->y(j+1))/
                (coordinates_m->y(j+2) - coordinates_m->y(j+1));
    f_xy[0] = (f_x[1][0] - f_x[0][0])*dy + f_x[0][0];
    f_xy[1] = (f_x[1][1] - f_x[0][1])*dy + f_x[0][1];
    // interpolation in z
    Value[0]  = (f_xy[1] - f_xy[0])/
                (coordinates_m->z(k+2) - coordinates_m->z(k+1))*
                (Point[2] - coordinates_m->z(k+1)) + f_xy[0];
}
}

