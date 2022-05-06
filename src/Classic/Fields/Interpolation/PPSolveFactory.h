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

#include "Fields/Interpolation/Mesh.h"
#include "Fields/Interpolation/PolynomialCoefficient.h"
#include "Fields/Interpolation/SquarePolynomialVector.h"

namespace interpolation {

class PolynomialPatch;

/** \class PPSolveFactory
 *  \brief PPSolveFactory solves the system of linear equations to interpolate
 *  from a grid of points using higher order polynomials, creating a
 *  PolynomialPatch object. The order of the polynomial is user-defined.
 *
 *  The PPSolveFactory can be used to match nearby points for smoothing or for
 *  straight fitting.
 *
 *  PPSolveFactory generates polynomials on a grid that is at the centre point
 *  of each of the mesh squares; PPSolveFactory picks nearby points in a
 *  square around the PPSolveFactory, picked by getNearbyPointsSquares. Two sets
 *  of points are chosen; the nearest points are used for fitting; smoothing is
 *  performed by choosing the points on the edge of the grid, and smoothing to
 *  derivatives on these points.
 *
 *  Because we work in a square grid, we have matching polynomials and
 *  derivatives; so e.g. 2D points are (0, 0), (1, 0), (0, 1), (1, 1) so
 *  corresponding polynomial coefficients that we solve for are
 *  x^0 y^0, x^1 y^0, x^0 y^1, x^1 y^1; note that hence we define "first order"
 *  to be all terms with no single x/y vector power > 1; i.e. "first order" in
 *  x OR first order in y, etc. We allow products so long as the maximum power
 *  in any one dimension <= 1
 *
 *  The PPSolveFactory sits on top of SolveFactory; PPSolveFactory has the job
 *  of picking values and derivatives for fitting, SolveFactory then does the
 *  actual solve for each individual polynomial.
 */
class PPSolveFactory {
  public:
    /** Constructor
     *
     *  \param points Set of points on which values are stored. Must be a
     *                rectangular grid. PPSolveFactory takes ownership of
     *                points (will delete on exit). The solve will return a
     *                wrong answer if the grid does not have regular spacing.
     *  \param values Set of values to which we fit. Must be one value per mesh
     *                point and each value must have the same size.
     *  \param polyPatchOrder The order of the fitted part of the polynomial.
     *  \param smoothingOrder The total order of the fitted and the smoothed
     *                 part of the polynomial. Smoothing order should always be
     *                 >= polyPatchOrder. So for example if polyPatchOrder
     *                 is 2 and smoothing order is 2, we get a 2nd order
     *                 function with no derivative matching. If polyPatchOrder
     *                 is 2 and smoothing order is 3, we get a 3rd order
     *                 function with derivative matching in first order at two
     *                 mesh points from the centre.
     */
    PPSolveFactory(Mesh* points,
                   std::vector<std::vector<double> > values,
                   int polyPatchOrder,
                   int smoothingOrder);

    /** Destructor deletes points_ */
    ~PPSolveFactory() {}

    /** Solve the system of equations to generate a PolynomialPatch object.
     *
     *  Returns a PolynomialPatch object, caller owns the returned memory.
     */
    PolynomialPatch* solve();

    /** Get nearby points in a square pattern
     *
     *  Get nearby points at least poly_order_lower distance away and less than
     *  poly_order_upper distance away; e.g. getNearbyPointsSquares(2, 4, 6)
     *  will get nearby points on a 2D grid at least 4 grid points and at most
     *  6 grid points, like:
     *
     *  x x x o o x ...
     *  x x x o o x
     *  x x x o o x
     *  o o o o o x
     *  o o o o o x
     *  x x x x x x
     *  ...
     */
    static std::vector<std::vector<int> > getNearbyPointsSquares(
                                        int pointDim,
                                        int polyOrderLower,
                                        int polyOrderUpper);

  private:
    void getPoints();
    void getValues(Mesh::Iterator it);
    void getDerivPoints();
    void getDerivs(const Mesh::Iterator& it);

    // nothing calls this method but I don't quite feel brave enough to remove
    // it...
    std::vector<double> outOfBoundsPosition(Mesh::Iterator outOfBoundsIt);
    static void nearbyPointsRecursive(
                                std::vector<int> check,
                                size_t checkIndex,
                                size_t polyPower,
                                std::vector<std::vector<int> >& nearbyPoints);

    PolynomialCoefficient getDeltaIterator(Mesh::Iterator it1,
                                           Mesh::Iterator it2, int valueIndex);

    int polyPatchOrder_m;
    int smoothingOrder_m;
    int polyDim_m;
    int valueDim_m;
    Mesh* polyMesh_m;
    Mesh* points_m;
    std::vector<std::vector<double> > values_m;
    std::vector<SquarePolynomialVector*> polynomials_m;

    std::vector< std::vector<double> > thisPoints_m;
    std::vector< std::vector<double> > thisValues_m;
    std::vector< std::vector<double> > derivPoints_m;
    std::vector< std::vector<double> > derivValues_m;
    std::vector< std::vector<int> > derivOrigins_m;
    std::vector< std::vector<int> > derivIndices_m;
    std::vector< MVector<double> > derivPolyVec_m;

    std::vector<std::vector<std::vector<int> > > edgePoints_m;
    std::vector< std::vector<int> > smoothingPoints_m;

};

}

