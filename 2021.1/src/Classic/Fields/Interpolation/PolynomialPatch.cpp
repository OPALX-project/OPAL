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

#include "Utilities/GeneralClassicException.h"

#include "Fields/Interpolation/SolveFactory.h"
#include "Fields/Interpolation/PolynomialPatch.h"

namespace interpolation {
PolynomialPatch::PolynomialPatch(Mesh* grid_points,
                    Mesh* validity_region,
                    std::vector<SquarePolynomialVector*> polynomials) {
    grid_points_ = grid_points;
    validity_region_ = validity_region;
    points_ = polynomials;
    // validate input data
    if (grid_points_ == NULL)
        throw GeneralClassicException(
                "PolynomialPatch::PolynomialPatch",
                "PolynomialPatch grid_points_ was NULL");
    if (validity_region_ == NULL)
        throw GeneralClassicException(
                "PolynomialPatch::PolynomialPatch",
                "PolynomialPatch validity_region_ was NULL");
    if (points_.size() == 0) {
        throw GeneralClassicException(
            "PolynomialPatch::PolynomialPatch",
            "Could not make PolynomialPatch with 0 length polynomials vector");
    }
    Mesh::Iterator end = grid_points_->end();
    if (int(points_.size()) != end.toInteger()) {
        throw GeneralClassicException(
            "PolynomialPatch::PolynomialPatch",
            "Could not make PolynomialPatch with bad length polynomials vector"
            );
    }
    point_dimension_ = grid_points_->getPositionDimension();
    if (validity_region_->getPositionDimension() != int(point_dimension_)) {
        throw GeneralClassicException(
                          "PolynomialPatch::PolynomialPatch",
                          "PolynomialPatch validity_region_ has bad dimension");
    }
    value_dimension_ = points_[0]->ValueDimension();
    for (size_t i = 0; i < points_.size(); ++i) {
        if (points_[i] == NULL)
            throw GeneralClassicException(
                              "PolynomialPatch::PolynomialPatch",
                              "PolynomialPatch points_ element was NULL");
        if (points_[i]->PointDimension() != point_dimension_) {
            throw GeneralClassicException(
                "PolynomialPatch::PolynomialPatch",
                "Polynomial with mismatched PointDimension in PolynomialPatch");
        }
        if (points_[i]->ValueDimension() != value_dimension_)
            throw GeneralClassicException(
                "PolynomialPatch::PolynomialPatch",
                "Polynomial with mismatched ValueDimension in PolynomialPatch");
    }
}

VectorMap* PolynomialPatch::clone() const {
    Mesh* new_mesh = NULL;
    if (grid_points_ != NULL) {
        new_mesh = grid_points_->clone();
    }
    Mesh* new_validity = NULL;
    if (validity_region_ != NULL) {
        new_validity = validity_region_->clone();
    }
    std::vector<SquarePolynomialVector*> new_points(points_.size());
    for (size_t i = 0; i < points_.size(); ++i) {
        if (points_[i] == NULL) {
            new_points[i] = NULL;
        } else {
            new_points[i] = new SquarePolynomialVector(*points_[i]);
        }
    }
    return new PolynomialPatch(new_mesh, new_validity, new_points);
}

PolynomialPatch::PolynomialPatch()
  : validity_region_(NULL), grid_points_(NULL), points_(), point_dimension_(0),
    value_dimension_(0) {
}

PolynomialPatch::~PolynomialPatch() {
    delete grid_points_;
    delete validity_region_; 
    for (size_t i = 0; i < points_.size(); ++i)
        delete points_[i];
}

void PolynomialPatch::function(const double* point, double* value) const {
    Mesh::Iterator nearest = grid_points_->getNearest(point);
    std::vector<double> point_temp(point_dimension_);
    std::vector<double> nearest_pos = nearest.getPosition();
    for (size_t i = 0; i < point_dimension_; ++i)
        point_temp[i] = point[i] - nearest_pos[i];
    int points_index = nearest.toInteger();
    points_[points_index]->F(&point_temp[0], value);
}

SquarePolynomialVector* PolynomialPatch::getPolynomialVector(const double* point) const {
    Mesh::Iterator nearest = grid_points_->getNearest(point);
    int points_index = nearest.toInteger();
    return points_[points_index];
}
}

