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

#include <vector>
#include <iostream>
#include <string>
#include <sstream>

#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_gamma.h>

#include "Utility/Inform.h" // ippl

#include "Utilities/GeneralClassicException.h"

#include "Fields/Interpolation/SolveFactory.h"
#include "Fields/Interpolation/PolynomialPatch.h"
#include "Fields/Interpolation/PPSolveFactory.h"

extern Inform* gmsg;

namespace interpolation {

typedef PolynomialCoefficient Coeff;

namespace {
    std::vector<double> getDeltaPos(Mesh* mesh) {
        Mesh::Iterator it = mesh->begin();
        std::vector<double> pos1 = it.getPosition();
        for (size_t i = 0; i < it.getState().size(); ++i)
            it[i]++;
        std::vector<double> pos2 = it.getPosition();
        for (size_t i = 0; i < pos2.size(); ++i)
            pos2[i] -= pos1[i];
        return pos2;
    }
}

PPSolveFactory::PPSolveFactory(Mesh* points,
                               std::vector<std::vector<double> > values,
                               int polyPatchOrder,
                               int smoothingOrder) 
  : polyPatchOrder_m(polyPatchOrder),
    smoothingOrder_m(smoothingOrder),
    polyDim_m(0),
    polyMesh_m(),
    points_m(points),
    values_m(values),
    polynomials_m(),
    thisPoints_m(),
    thisValues_m(),
    derivPoints_m(),
    derivValues_m(),
    derivIndices_m(),
    edgePoints_m(),
    smoothingPoints_m() {
    if (points == nullptr) {
        throw GeneralClassicException(
                          "PPSolveFactory::Solve",
                          "nullptr Mesh input to Solve");
    }
    if (points->end().toInteger() != int(values.size())) {
        std::stringstream ss;
        ss << "Mismatch between Mesh and values in Solve. Mesh size was "
           << points->end().toInteger() << " but field values size was "
           << values.size();
        throw GeneralClassicException(
                          "PPSolveFactory::Solve",
                          ss.str());
    }
    if (smoothingOrder < polyPatchOrder) {
        throw GeneralClassicException(
                "PPSolveFactory::Solve",
                "Polynomial order must be <= smoothing order in Solve");
    }
    polyMesh_m = points->dual();
    if (polyMesh_m == nullptr)
        throw GeneralClassicException(
                          "PPSolveFactory::Solve",
                          "Dual of Mesh was nullptr");
    polyDim_m = values[0].size();
    int posDim = points->getPositionDimension();

    edgePoints_m = std::vector< std::vector< std::vector<int> > >(posDim+1);
    for (int i = 1; i <= posDim; ++i)
        edgePoints_m[i] = getNearbyPointsSquares
                                      (i, 0, smoothingOrder-polyPatchOrder);
    smoothingPoints_m = getNearbyPointsSquares
                                (posDim, polyPatchOrder_m-1, polyPatchOrder_m);
}

std::vector<double> PPSolveFactory::outOfBoundsPosition(
                                              Mesh::Iterator outOfBoundsIt) {
    const Mesh* mesh = outOfBoundsIt.getMesh();
    std::vector<int> state = outOfBoundsIt.getState();
    std::vector<int> delta = std::vector<int>(state.size(), 2);
    Mesh::Iterator end = mesh->end()-1;
    std::vector<double> pos1 = mesh->begin().getPosition();
    std::vector<double> pos2 = Mesh::Iterator(delta, mesh).getPosition();
    for (size_t i = 0; i < state.size(); ++i)
        if (state[i] < 1) {
            state[i] = 1;
        } else if (state[i] > end[i]) { 
            state[i] = end[i];
        }
    std::vector<double> position = Mesh::Iterator(state, mesh).getPosition();
    state = outOfBoundsIt.getState();
    for (size_t i = 0; i < state.size(); ++i)
        if (state[i] < 1) {
            position[i] = (pos2[i] - pos1[i])*(state[i]-1) + position[i];
        } else if (state[i] > end[i]) { 
            position[i] = (pos2[i] - pos1[i])*(state[i]-end[i]) + position[i];
        }
    return position;
}

void PPSolveFactory::getPoints() {
    int posDim = points_m->getPositionDimension();
    std::vector< std::vector<int> > dataPoints =
                          getNearbyPointsSquares(posDim, -1, polyPatchOrder_m);
    thisPoints_m = std::vector< std::vector<double> >(dataPoints.size());
    std::vector<double> deltaPos = getDeltaPos(points_m);
    // now iterate over the index_by_power_list elements - offset; each of
    // these makes a point in the polynomial fit
    for (size_t i = 0; i < dataPoints.size(); ++i) {
        thisPoints_m[i] = std::vector<double>(posDim);
        for (int j = 0; j < posDim; ++j)
            thisPoints_m[i][j] = (0.5-dataPoints[i][j])*deltaPos[j];
     }
}

void PPSolveFactory::getValues(Mesh::Iterator it) {
    thisValues_m = std::vector< std::vector<double> >();
    int posDim = it.getState().size();
    std::vector< std::vector<int> > dataPoints =
                          getNearbyPointsSquares(posDim, -1, polyPatchOrder_m);
    size_t dataPointsSize = dataPoints.size();
    Mesh::Iterator end = points_m->end()-1;
    // now iterate over the indexByPowerList elements - offset; each of
    // these makes a point in the polynomial fit
    for (size_t i = 0; i < dataPointsSize; ++i) {
        std::vector<int> itState = it.getState();
        for (int j = 0; j < posDim; ++j)
            itState[j]++;
        Mesh::Iterator itCurrent = Mesh::Iterator(itState, points_m);
        bool outOfBounds = false; // element is off the edge of the mesh
        for (int j = 0; j < posDim; ++j) {
            itCurrent[j] -= dataPoints[i][j];
            outOfBounds = outOfBounds ||
                            itCurrent[j] < 1 ||
                            itCurrent[j] > end[j];
        }
        if (outOfBounds) { // if off the edge, then just constrain to zero
            thisValues_m.push_back(values_m[it.toInteger()]);
        } else { // else fit using values
            thisValues_m.push_back(values_m[itCurrent.toInteger()]);
        }
     }
}

void PPSolveFactory::getDerivPoints() {
    derivPoints_m = std::vector< std::vector<double> >();
    derivIndices_m = std::vector< std::vector<int> >();
    derivOrigins_m = std::vector< std::vector<int> >();
    derivPolyVec_m = std::vector<MVector<double> >();
    int posDim = points_m->getPositionDimension();
    // get the outer layer of points
    int deltaOrder = smoothingOrder_m - polyPatchOrder_m;
    if (deltaOrder <= 0)
        return;
    std::vector<double> deltaPos = getDeltaPos(points_m);
    // smoothingPoints_m is the last layer of points used for polynomial fit
    // walk around smoothingPoints_m finding the points used for smoothing
    for (size_t i = 0; i < smoothingPoints_m.size(); ++i) {
        // make a list of the axes that are on the edge of the space
        std::vector<int> equalAxes;
        for (size_t j = 0; j < smoothingPoints_m[i].size(); ++j)
            if (smoothingPoints_m[i][j] == polyPatchOrder_m)
                equalAxes.push_back(j);
        std::vector<std::vector<int> > edgePoints =
                                                edgePoints_m[equalAxes.size()];
        // note the first point, 0,0, is ignored
        for (size_t j = 0; j < edgePoints.size(); ++j) {
            derivIndices_m.push_back(std::vector<int>(posDim));
            derivOrigins_m.push_back(smoothingPoints_m[i]);
            derivPoints_m.push_back(std::vector<double>(posDim));

            for (size_t k = 0; k < smoothingPoints_m[i].size(); ++k) {
                derivPoints_m.back()[k] = (smoothingPoints_m[i][k]-0.5)*deltaPos[k];
            }
            for (size_t k = 0; k < edgePoints[j].size(); ++k) {
                derivIndices_m.back()[equalAxes[k]] = edgePoints[j][k];
            }
        }
    }

    int nPolyCoeffs = SquarePolynomialVector::NumberOfPolynomialCoefficients
                                                    (posDim, smoothingOrder_m);
    for (size_t i = 0; i < derivPoints_m.size(); ++i) {
        derivPolyVec_m.push_back(MVector<double>(nPolyCoeffs, 1.));
        // Product[x^{p_j-q_j}*p_j!/(p_j-q_j)!]
        for (int j = 0; j < nPolyCoeffs; ++j) {
            std::vector<int> pVec =
                               SquarePolynomialVector::IndexByPower(j, posDim);
            std::vector<int> qVec = derivIndices_m[i];
            for (int l = 0; l < posDim; ++l) {
                int pow = pVec[l]-qVec[l];
                if (pow < 0) {
                    derivPolyVec_m.back()(j+1) = 0.;
                    break;
                }
                derivPolyVec_m.back()(j+1) *= gsl_sf_fact(pVec[l])/gsl_sf_fact(pow)*gsl_sf_pow_int(-deltaPos[l]/2., pow);
            }
        }
    }
}

void PPSolveFactory::getDerivs(const Mesh::Iterator& it) {
    // calculate the derivatives near to polynomial at position indexed by it
    #ifdef DEBUG
    bool verbose = false;
    if (verbose) {
        std::cerr << "PPSolveFactory::GetDerivs at position ";
        for (auto p: it.getPosition()) {std::cerr << p << " ";}
        std::cerr << std::endl;
    }
    #endif
    // derivatives go in derivValues_m
    derivValues_m = std::vector< std::vector<double> >(derivPoints_m.size());

    std::vector<double> itPos = it.getPosition();
    int posDim = it.getState().size();

    // get the outer layer of points
    Mesh::Iterator last = it.getMesh()->end()-1;
    int deltaOrder = smoothingOrder_m - polyPatchOrder_m;
    if (deltaOrder <= 0) { // no smoothing needed
        return;
    }
    for (size_t i = 0; i < derivPoints_m.size(); ++i) {
        Mesh::Iterator nearest = it;
        bool inTheMesh = true;
        for (int j = 0; j < posDim && inTheMesh; ++j) {
            nearest[j] += derivOrigins_m[i][j];
            inTheMesh = nearest[j] <= last[j];
        }
        if (!inTheMesh) {
            derivValues_m[i] = std::vector<double>(polyDim_m, 0.);
        } else {
            MMatrix<double> coeffs =
                  polynomials_m[nearest.toInteger()]->GetCoefficientsAsMatrix();
            MVector<double> values = coeffs*derivPolyVec_m[i];
            derivValues_m[i] = std::vector<double>(polyDim_m);
            for(int j = 0; j < posDim; ++j) {
                derivValues_m[i][j] = values(j+1);
            }
        }
        #ifdef DEBUG
        if (verbose) {
            std::cerr << "Point ";
            for (auto p: derivPoints_m[i]) {std::cerr << p << " ";}
            std::cerr << " values ";
            for (auto v: derivValues_m[i]) {std::cerr << v << " ";}
            std::cerr << " derivIndices ";
            for (auto in: derivIndices_m[i]) {std::cerr << in << " ";}
            std::cerr << std::endl;
        }
        #endif
    }
}

PolynomialPatch* PPSolveFactory::solve() {
    Mesh::Iterator begin = points_m->begin();
    Mesh::Iterator end = points_m->end()-1;
    int meshSize = polyMesh_m->end().toInteger();
    polynomials_m = std::vector<SquarePolynomialVector*>(meshSize, nullptr);
    // get the list of points that are needed to make a given poly vector
    getPoints();
    getDerivPoints();
    SolveFactory solver(smoothingOrder_m,
                        polyMesh_m->getPositionDimension(),
                        polyDim_m,
                        thisPoints_m,
                        derivPoints_m,
                        derivIndices_m);
    int total = (polyMesh_m->end()-1).toInteger();
    double oldPercentage = 0.;
    for (Mesh::Iterator it = polyMesh_m->end()-1;
         it >= polyMesh_m->begin();
         --it) {
        double newPercentage = (total-it.toInteger())/double(total)*100.;
        if (newPercentage - oldPercentage > 10.) {
            *gmsg << "    Done " << newPercentage << " %" << endl;
            oldPercentage = newPercentage;
        }
        // find the set of values and derivatives for this point in the mesh
        getValues(it);
        getDerivs(it);
        // The polynomial is found using simultaneous equation solve
        polynomials_m[it.toInteger()] = solver.PolynomialSolve
                                                  (thisValues_m, derivValues_m);
    }
    return new PolynomialPatch(polyMesh_m, points_m, polynomials_m);
}

void PPSolveFactory::nearbyPointsRecursive(
                               std::vector<int> check,
                               size_t checkIndex,
                               size_t polyPower,
                               std::vector<std::vector<int> >& nearbyPoints) {
    check[checkIndex] = polyPower;
    nearbyPoints.push_back(check);
    if (checkIndex+1 == check.size())
        return;
    for (int i = 1; i < int(polyPower); ++i)
        nearbyPointsRecursive(check, checkIndex+1, i, nearbyPoints);
    for (int i = 0; i <= int(polyPower); ++i) {
        check[checkIndex] = i;
        nearbyPointsRecursive(check, checkIndex+1, polyPower, nearbyPoints);
    }
}

std::vector<std::vector<int> > PPSolveFactory::getNearbyPointsSquares(
                                                      int pointDim,
                                                      int polyOrderLower,
                                                      int polyOrderUpper) {
    if (pointDim < 1)
        throw(GeneralClassicException("PPSolveFactory::getNearbyPointsSquares",
                                      "Point dimension must be > 0"));
    if (polyOrderLower > polyOrderUpper)
        throw(GeneralClassicException("PPSolveFactory::getNearbyPointsSquares",
                   "Polynomial lower bound must be <= polynomial upper bound"));
    // polyOrder -1 (or less) means no terms
    if (polyOrderUpper == polyOrderLower || polyOrderUpper < 0)
        return std::vector<std::vector<int> >();
    // polyOrder 0 means constant term
    std::vector<std::vector<int> > nearbyPoints(1, std::vector<int>(pointDim, 0));
    // polyOrder > 0 means corresponding poly terms (1 is linear, 2 is quadratic...)
    for (size_t thisPolyOrder = 0;
         thisPolyOrder < size_t(polyOrderUpper);
         ++thisPolyOrder) 
        nearbyPointsRecursive(nearbyPoints[0], 0, thisPolyOrder+1, nearbyPoints);
    if (polyOrderLower < 0)
        return nearbyPoints; // no terms in lower bound
    int nPointsLower = 1;
    for (int dim = 0; dim < pointDim; ++dim)
        nPointsLower *= polyOrderLower+1;
    nearbyPoints.erase(nearbyPoints.begin(), nearbyPoints.begin()+nPointsLower);
    return nearbyPoints;
}
}
