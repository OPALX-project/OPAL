/*
 *  Copyright (c) 2012-2014, Chris Rogers
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

#ifndef RING_SECTION_H
#define RING_SECTION_H

#include <vector>

#include "AbsBeamline/Component.h"

/** \class[RingSection]
 *
 *  \brief Component placement handler in ring geometry
 *
 *  RingSection handles placement of a component when it is placed in a ring
 *  geometry. Here the primary index for section placement is azimuthal angle.
 *  RingSection assumes that ring objects occupy a space defined by a start
 *  plane and another end plane.
 *
 *  All vectors should be in cartesian coordinates (x,y,z); with z being the
 *  axis of the ring.
 *
 *  \param component_m component that the OpalSection wraps - this is a borrowed
 *                   reference (RingSection does not own the memory)
 *  \param startPosition_m position of the centre of the start face in
 *                   cylindrical polar coordinates
 *  \param startOrientation_m vector normal to the start face pointing towards
 *                   the field map
 *  \param endPosition_m position of the centre of the end face in cylindrical
 *                   polar coordinates
 *  \param endOrientation_m vector normal to the end face pointing away from the
 *                   field map
 *  \param componentPosition_m field map position relative to component start
 *  \param componentRotation_m field map rotation R
 *
 *  So field maps are calculated in local coordinate system
 *          U_local = V+R*U_global
 *  where R is rotation matrix and V is componentPosition_m.
 *
 *  Return field values are returned like B_global = R^{-1}*B_local
 */

class RingSection {
  public:
    /** Construct a ring section - positions, orientations etc default to 0.
     */
    RingSection();

    /** Copy constructor; deepcopies the Component (and copies everything else)*/
    RingSection(const RingSection& sec);

    RingSection& operator=(const RingSection& sec);

    /** Destructor - does nothing*/
    ~RingSection();

    /** Return true if pos is on or past start plane
     *
     *  \param   pos position to test
     *  \returns true if phi(pos) >= phi(start plane), after taking account of
     *           the position and rotation of the start face
     */
    bool isOnOrPastStartPlane(const Vector_t& pos) const;

    /** Return true if pos is past end plane
     *
     *  \returns true if phi > phi(end plane), after taking account of the
     *           position and rotation of the start face
     */
    bool isPastEndPlane(const Vector_t& pos) const;

    /** Return field value in global coordinate system
     *
     *  \param pos Position in global Cartesian coordinates
     *  \param t time in lab frame
     *  \param centroid Not sure what this is
     *  \param E Vector to be filled with electric field values; will always
     *         overwrite with 0. before filling field
     *  \param B Vector to be filled with magnetic field values; will always
     *         overwrite with 0. before filling field
     *  \returns true if pos is outside of the field bounding box, else false
     *
     *  Note; does not check for component == NULL; caller must assign
     *  component_m (using setComponent) before calling this function.
     */
    bool getFieldValue(const Vector_t& pos, const Vector_t& centroid,
                       const double& t, Vector_t& E, Vector_t& B) const;

    /** Get the "Virtual" bounding box for the RingSection
     *
     *  Defined by the point one radial distance along each face of the element
     *  either towards the centre of the ring or towards the outside of the
     *  ring
     *
     *  Return order is; start outside, start inside, end outside, end inside
     *
     *  An improvement would be to put the crossing point if the two faces cross
     *  before they reach the radial apertures described above.
     */
    std::vector<Vector_t> getVirtualBoundingBox() const;

    /** Return true if the phi range overlaps bounding box elements */
    bool doesOverlap(double phiStart, double phiEnd) const;

    /** Set the component wrapped by RingSection
     *
     *  This borrows the Component* pointer (caller is responsible for cleanup)
     */
    inline void setComponent(Component* component) {component_m = component;}

    /** Get the component wrapped by RingSection
     *
     *  Component* is not owned by caller or RingSection
     */
    inline Component* getComponent() const {return component_m;}

    /** Set a position on the plane of the section start */
    inline void setStartPosition(Vector_t pos) {startPosition_m = pos;}

    /** Get a position on the plane of the section start */
    inline Vector_t getStartPosition() const {return startPosition_m;}

    /** Set the normal vector to the section start plane */
    inline void setStartNormal(Vector_t orientation);

    /** Get the normal vector to the section start plane */
    inline Vector_t getStartNormal() const {return startOrientation_m;}

    /** Set a position on the section end plane */
    inline void setEndPosition(Vector_t pos) {endPosition_m = pos;}

    /** Get a position on the section end plane */
    inline Vector_t getEndPosition() const {return endPosition_m;}

    /** Set the normal vector to the section end plane  */
    inline void setEndNormal(Vector_t orientation);

    /** Get the normal vector to the section end plane */
    inline Vector_t getEndNormal() const {return endOrientation_m;}

    /** Set the displacement for the component relative to the section start */
    inline void setComponentPosition(Vector_t position) {componentPosition_m = position;}

    /** Get the displacement for the component relative to the section start */
    inline Vector_t getComponentPosition() const {return componentPosition_m;}

    /** Set the rotation for the component relative to the section start */
    inline void setComponentOrientation(Vector_t orientation);

    /** Get the rotation for the component relative to the section start */
    inline Vector_t getComponentOrientation() const {return componentOrientation_m;}
  private:
    void rotate(Vector_t& vector) const;
    void rotate_back(Vector_t& vector) const;
    inline Vector_t& normalise(Vector_t& vector) const;
    inline void rotateToTCoordinates(Vector_t& vec) const;
    inline void rotateToCyclCoordinates(Vector_t& vec) const;

    Component* component_m;

    Vector_t componentPosition_m;
    Vector_t componentOrientation_m;

    Vector_t startPosition_m;
    Vector_t startOrientation_m;

    Vector_t endPosition_m;
    Vector_t endOrientation_m;

    void updateComponentOrientation();
    double sin2_m;
    double cos2_m;
};

typedef std::vector<RingSection*> RingSectionList;

inline void RingSection::setComponentOrientation(Vector_t orientation) {
    componentOrientation_m = orientation;
    updateComponentOrientation();
}

inline void RingSection::setStartNormal(Vector_t orientation) {
    startOrientation_m = normalise(orientation);
}

inline void RingSection::setEndNormal(Vector_t orientation) {
    endOrientation_m = normalise(orientation);
}

inline Vector_t& RingSection::normalise(Vector_t& orientation) const {
    double magnitude = sqrt(orientation(0)*orientation(0)+
                            orientation(1)*orientation(1)+
                            orientation(2)*orientation(2));
    if (magnitude > 0.)
        orientation /= magnitude;
    return orientation;
}


void RingSection::rotateToTCoordinates(Vector_t& vec) const {
  vec = Vector_t(vec(1), vec(2), vec(0));
}

void RingSection::rotateToCyclCoordinates(Vector_t& vec) const {
  vec = Vector_t(vec(2), vec(0), vec(1));
}

#endif //RING_SECTION_H
