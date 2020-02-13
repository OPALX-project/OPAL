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

#include "Fields/SectorField.h"

#include <math.h>

#include <limits>
#include <vector>
#include <algorithm>

#include "Utilities/LogicalError.h"

// note I don't use exactly max() as the bounding box because I don't want to
// accidentally wrap around due to some unforeseen floating point precision
// issue
SectorField::SectorField(const std::string& /*file_name*/)
    : bbMin_m(3, -std::numeric_limits<double>::max()/10.),
      bbMax_m(3, std::numeric_limits<double>::max()/10.),
      polarBBMin_m(3, 0),  polarBBMax_m(3, 0) {
    polarBBMin_m[1] = bbMin_m[1];
    polarBBMax_m[1] = bbMax_m[1];
    polarBBMax_m[0] = bbMax_m[0];
    polarBBMin_m[2] = -2.*M_PI;
    polarBBMax_m[2] = 2.*M_PI;
}

SectorField::~SectorField() {}

void SectorField::convertToPolar(double* position) {
    double x = ::sqrt(position[0]*position[0]+position[2]*position[2]);
    double z = ::atan2(position[2], position[0]);
    position[0] = x;
    position[2] = z;
}

void SectorField::convertToPolar(const double* position, double* value) {
    double x = +value[0]*::cos(position[2])
               +value[2]*::sin(position[2]);
    double z = +value[2]*::cos(position[2])
               -value[0]*::sin(position[2]);
    value[0] = x;
    value[2] = z;
}

void SectorField::convertToCartesian(double* position) {
    double x = position[0]*::cos(position[2]);  // r cos(phi)
    double z = position[0]*::sin(position[2]);  // r sin(phi)
    position[0] = x;
    position[2] = z;
}

void SectorField::convertToCartesian(const double* position, double* value) {
    double x = +value[0]*::cos(position[2])
               -value[2]*::sin(position[2]);
    double z = +value[2]*::cos(position[2])
               +value[0]*::sin(position[2]);
    value[0] = x;
    value[2] = z;
}

void SectorField::setPolarBoundingBox
                        (double bbMinR, double bbMinY, double bbMinPhi,
                         double bbMaxR, double bbMaxY, double bbMaxPhi,
                         double bbTolR, double bbTolY, double bbTolPhi) {
    if (bbMinR > bbMaxR) {
        throw(LogicalError(
               "SectorField::SetPolarBoundingBox",
               "Bounding box minimum radius was greater than maximum radius"));
    }
    if (bbMinR < 0.) {
        throw(LogicalError(
               "SectorField::SetPolarBoundingBox",
               "Bounding box radius must be positive"
               ));
    }
    if (bbMinY > bbMaxY) {
        throw(LogicalError(
               "SectorField::SetPolarBoundingBox",
               "Bounding box minimum y was greater than maximum y"
               ));
    }
    if (bbMinPhi > bbMaxPhi) {
        throw(LogicalError(
               "SectorField::SetPolarBoundingBox",
               "Bounding box minimum angle was greater than maximum angle"
               ));
    }
    if (bbMinPhi < -2.*M_PI || bbMinPhi > 2.*M_PI ||
        bbMaxPhi < -2.*M_PI || bbMaxPhi > 2.*M_PI) {
        throw(LogicalError(
               "SectorField::SetPolarBoundingBox",
               "Bounding box angles must be in range -2*M_PI < phi < 2*M_PI"
               ));
    }

    polarBBMin_m[0] = bbMinR-bbTolR;
    polarBBMin_m[1] = bbMinY-bbTolY;
    polarBBMin_m[2] = bbMinPhi-bbTolPhi;

    polarBBMax_m[0] = bbMaxR+bbTolR;
    polarBBMax_m[1] = bbMaxY+bbTolY;
    polarBBMax_m[2] = bbMaxPhi+bbTolPhi;

    // bounding box from corner coordinates
    std::vector< std::vector<double> > corner_coords(
                                getCorners(bbMinR, bbMinPhi, bbMaxR, bbMaxPhi));
    bbMin_m[0] =
            *std::min_element(corner_coords[0].begin(), corner_coords[0].end());
    bbMax_m[0] =
            *std::max_element(corner_coords[0].begin(), corner_coords[0].end());
    bbMin_m[1] = bbMinY;
    bbMax_m[1] = bbMaxY;
    bbMin_m[2] =
            *std::min_element(corner_coords[1].begin(), corner_coords[1].end());
    bbMax_m[2] =
            *std::max_element(corner_coords[1].begin(), corner_coords[1].end());

    // if the magnet crosses an axis, then the corners are no longer at the max
    // extent
    if ( (bbMaxPhi > 0.5*M_PI && bbMinPhi < 0.5*M_PI) ||
         (bbMaxPhi > -1.5*M_PI && bbMinPhi < -1.5*M_PI) ) {
        bbMax_m[2] = bbMaxR;
    }
    if ((bbMaxPhi > M_PI && bbMinPhi < M_PI) ||
        (bbMaxPhi > -M_PI && bbMinPhi < -M_PI)) {
        bbMin_m[0] = -bbMaxR;
    }
    if ((bbMaxPhi > 1.5*M_PI && bbMinPhi < 1.5*M_PI) ||
        (bbMaxPhi > -0.5*M_PI && bbMinPhi < -0.5*M_PI)) {
        bbMin_m[2] = -bbMaxR;
    }
    if ((bbMaxPhi > 0.*M_PI && bbMinPhi < 0.*M_PI)) {
        bbMax_m[0] = bbMaxR;
    }
}

std::vector< std::vector<double> > SectorField::getCorners
              (double bbMinR, double bbMinPhi, double bbMaxR, double bbMaxPhi) {
    std::vector< std::vector<double> > corner_coords(2);
    corner_coords[0] = std::vector<double>(4);
    corner_coords[1] = std::vector<double>(4);
    // corners in polar coordinates
    double corner_0[3] = {bbMinR, 0., bbMinPhi};
    double corner_1[3] = {bbMinR, 0., bbMaxPhi};
    double corner_2[3] = {bbMaxR, 0., bbMaxPhi};
    double corner_3[3] = {bbMaxR, 0., bbMinPhi};
    convertToCartesian(corner_0);
    convertToCartesian(corner_1);
    convertToCartesian(corner_2);
    convertToCartesian(corner_3);
    // corners in rectangular coordinates (ignore y)
    corner_coords[0][0] = corner_0[0];
    corner_coords[0][1] = corner_1[0];
    corner_coords[0][2] = corner_2[0];
    corner_coords[0][3] = corner_3[0];
    corner_coords[1][0] = corner_0[2];
    corner_coords[1][1] = corner_1[2];
    corner_coords[1][2] = corner_2[2];
    corner_coords[1][3] = corner_3[2];
    return corner_coords;
}


std::vector<double> SectorField::getPolarBoundingBoxMin() const {
    return polarBBMin_m;
}

std::vector<double> SectorField::getPolarBoundingBoxMax() const {
    return polarBBMax_m;
}

void SectorField::getFieldDimensions(double &zBegin, double &zEnd,
                                     double &rBegin, double &rEnd) const {
    zBegin = polarBBMin_m[2]*(polarBBMin_m[0]+polarBBMax_m[0])/2.;
    zEnd = polarBBMax_m[2]*(polarBBMin_m[0]+polarBBMax_m[0])/2.;
    rBegin = polarBBMin_m[0];
    rEnd = polarBBMin_m[2];
}

void SectorField::getFieldDimensions(double &xIni, double &xFinal,
                                     double &yIni, double &yFinal,
                                     double &zIni, double &zFinal) const {
    xIni = bbMin_m[0];
    yIni = bbMin_m[1];
    zIni = bbMin_m[2];
    xFinal = bbMax_m[0];
    yFinal = bbMax_m[1];
    zFinal = bbMax_m[2];
}
