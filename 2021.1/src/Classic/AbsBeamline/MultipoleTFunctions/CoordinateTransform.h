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

#ifndef COORDINATE_TRANSFORM_H
#define COORDINATE_TRANSFORM_H

/** ---------------------------------------------------------------------
  *
  * CoordinateTransform transforms between the coordinate system centred \n 
  * at the entrance and the local Frenet-Serret coordinate system. Assuming \n
  * CoordinateTransform takes in the fringe field parameters and calculates \n
  * the coordinate transformation. First the coordinates are transformed to\n
  * a coordinate system centred in the middle of the magnet. From here \n
  * the Frenet-Serret coordinates are found using numerical integration.
  *
  * ---------------------------------------------------------------------
  *
  * Class category: AbsBeamline \n
  * Author: Martin Duy Tat\n
  *
  * ---------------------------------------------------------------------
  */

#include <vector>

namespace coordinatetransform {

struct myParams {
    double s_0;
    double lambdaleft;
    double lambdaright;
    double rho;
};

class CoordinateTransform {
public:
    /** Default constructor, transforms everything to the origin */
    CoordinateTransform() = delete;
    /** Constructor, calculates coordinate transformation from lab coordinates
     *  to Frenet-Serret coordinates, given fringe field parameters
     *  \param xlab -> x-coordinate in lab frame
     *  \param ylab -> y-coordinate in lab frame
     *  \param zlab -> z-coordinate in lab frame
     *  \param s_0 -> Centre field length
     *  \param lambdaleft -> Left end field length
     *  \param lambdaright -> Right end field length
     *  \param rho -> Centre radius of curvature
     */
    CoordinateTransform(const double &xlab,
                        const double &ylab,
                        const double &zlab,
                        const double &s_0,
                        const double &lambdaleft,
                        const double &lambdaright,
                        const double &rho);
    /** Copy constructor */
    CoordinateTransform(const CoordinateTransform &transform);
    /** Destructor, does nothing */
    ~CoordinateTransform();
    /** Assigment operator */
    CoordinateTransform& operator= (const CoordinateTransform &transform);
    /** Returns a list of transformed coordinates */
    std::vector<double> getTransformation() const;
    /** Calculates reference trajectory by integrating the unit tangent vector
     *  \param s -> s-coordinate in local Frenet-Serret coordinates
     */
    std::vector<double> calcReferenceTrajectory(const double &s) const;
    /** Returns unit tangent vector
     *  \param s -> s-coordinate in local Frenet-Serret coordinates
     */
    std::vector<double> getUnitTangentVector(const double &s) const;
private:
    /** Calculates the coordinate s
     *  \param xlab -> x-coordinate in lab frame
     *  \param ylab -> y-coordinate in lab frame
     */
    void calcSCoordinate(const double &xlab, const double &ylab);
    /** Calculates the coordinate x, coordinate s must be calculated first! 
     *  \param xlab -> x-coordinate in lab frame
     *  \param ylab -> y-coordinate in lab frame
     */
    void calcXCoordinate(const double &xlab, const double &ylab);
    /** Transforms from coordinate system centred in the middle of the magnet
     *  to the coordinate system placed at the entrance
     *  \param coordinates -> Coordinates in coordinate system centred
     *  in the middle of the magnet
     *  \param boundingBoxLength -> Length along the magnet from the magnet
     *  entrance to the middle of the magnet
     */
    void transformFromEntranceCoordinates(std::vector<double> &coordinates,
                                          const double &boundingBoxLength);
    double s_0_m;
    double lambdaleft_m;
    double lambdaright_m;
    double rho_m;
    double x_m;
    double z_m;
    double s_m;
    static const double error;
    static const int workspaceSize;
    static const int algorithm;
};

/** Internal function used for GSL numerical integration */
double getUnitTangentVectorX(double s, void *p);
/** Internal function used for GSL numerical integration */
double getUnitTangentVectorY(double s, void *p);

}
#endif
