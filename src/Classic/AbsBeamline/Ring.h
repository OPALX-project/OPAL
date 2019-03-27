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

#ifndef RING_H
#define RING_H

#include <string>

#include "AbsBeamline/Component.h"

#include "BeamlineGeometry/PlanarArcGeometry.h"

#include "Utilities/RingSection.h"
#include "Utilities/GeneralClassicException.h"

class LossDataSink;
template <class T, unsigned Dim>
class PartBunchBase;
class FieldMap;

/** \class[Ring]
 *
 *  \brief Ring describes a ring type geometry for tracking
 *
 *  Ring describes a ring type geometry for tracking. Ring provides the
 *  necessary interfaces for e.g. OPAL-CYCL to track through the ring geometry,
 *  while enabling the user to add arbitrary field elements in a closed
 *  geometry.
 *
 *  Ring uses similar routines to OPAL-T OpalBeamline class to set up
 *  geometry; note that as OPAL-CYCL places the beam in an x-y geometry we place
 *  in an x-y geometry for Ring (i.e. the axis of the ring is by default the
 *  z-axis). So far only placement of elements in the midplane is supported. It
 *  is not possible to give vertical displacements or rotations, or add elements
 *  that might create vertical displacements and rotations.
 *
 *  Also aim to maintain backwards compatibility with Cyclotron (i.e. use
 *  ParallelCyclotronTracker)
 */
class Ring : public Component {
  public:
    /** Constructor
     *
     *  \param ring Name of the ring as defined in the input file
     */
    Ring(std::string ring);

    /** Copy constructor
     *
     *  Can't copy LossDataSink so throw exception if this is set
     */
    Ring(const Ring& ring);

    /** Destructor - deletes lossDS_m if not NULL */
    virtual ~Ring();

    /** Overwrite data in vector E and B with electric and magnetic field
     *
     *  @param i index of item in RefPartBunch_m - particle bunch
     *  @param t time
     *  @param E array where electric field vector is stored - any
     *         existing data is overwritten
     *  @param B array where magnetic field vector is stored - any
     *         existing data is overwritten
     *
     *  @returns false if particle is outside of the field map apertures, else
     *  true. If particle is off the field maps, then set flag on the particle
     *  "Bin" data to -1 and store in LossDataSink
     */
    virtual bool apply(const size_t &id, const double &t, Vector_t &E,
                       Vector_t &B) override;

    /** Overwrite data in vector E and B with electromagnetic field at point R
     *
     *  @param R 3 vector position at which the field is found in Cartesian
     *         coordinates (i.e. x, y, z with z=vertical)
     *  @param P 3 vector momentum
     *  @param centroid unknown, but not used - bunch mean maybe?
     *  @param t time
     *  @param E vector where electric field vector will be stored - any
     *         existing data is overwritten
     *  @param B vector where magnetic field vector will be stored - any
     *         existing data is overwritten
     *
     *  @returns false if particle is outside of the ring aperture, else true.
     *  If particle is off the field maps, then set flag on the particle
     *  "Bin" data to -1
     */
    virtual bool apply(const Vector_t &R, const Vector_t &P,
                       const double &t, Vector_t &E, Vector_t &B) override;

    /** Initialise the Ring
     *
     *  @param bunch the particle bunch. Ring borrows this pointer (caller
     *         owns memory)
     *  @param startField - not used
     *  @param endField - not used
     *  @param scaleFactor - not used
     */
    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField,
                            double &endField) override;

    /** Initialise the Ring - set the bunch and allocate a new LossDataSink
     *
     *  @param bunch the particle bunch. Ring borrows this pointer (caller
     *         owns memory)
     */
    virtual void initialise(PartBunchBase<double, 3> *bunch);

    /** Clean up the Ring
     *
     *  Ring relinquishes RefPartBunchBase<double, 3> pointer and deletes LossDataSink
     */
    virtual void finalise() override;

    /** Returns true - Ring is assumed to bend particles, being a ring */
    virtual bool bends() const override {return true;}

    /** Accept the BeamlineVisitor
     *
     *  Just calls visitRing function on the visitor. I guess the point of
     *  this function is that it enables us to store a pointer to the visitor
     *  object or something
     */
    virtual void accept(BeamlineVisitor& visitor) const override;

    /** Not implemented - always throws an exception */
    virtual void getDimensions(double &zBegin, double &zEnd) const override;

    /** Inherited copy constructor */
    virtual ElementBase* clone() const override {return new Ring(*this);}

    /** Add element to the ring
     *
     *  Add element to the ring. Elements are assumed to occupy a region of
     *  space defined by a (flat) plane at the start and a plane at the end,
     *  both infinite in extent. The position and rotation of these planes are
     *  defined according to the Component geometry given by
     *  element.getGeometry().
     *
     *  Caller owns memory allocated to element - Ring makes a copy.
     *
     *  Throws an exception if the geometry would bend the element out of the
     *  midplane (elements out of midplane are not yet supported). Note that
     *  BeamlineGeometry considers midplane to be x-z, whereas Ring
     *  considers midplane to be x-y; we apply a rotation during set up.
     *
     *  The element is assumed to extend not beyond the element geometry;
     *  Ring applies a bounding box based on the element geometry, if there
     *  are field maps expanding outside this region they may get cut.
     */
    void appendElement(const Component &element);

    /** Not implemented, throws an exception */
    virtual EMField &getField() override {throw GeneralClassicException("Ring::getField", "Not implemented");}

    /** Not implemented, throws an exception */
    virtual const EMField &getField() const override {throw GeneralClassicException("Ring::getField", "Not implemented");}

    /** Not implemented */
    virtual PlanarArcGeometry &getGeometry() override {return planarArcGeometry_m;}

    /** Not implemented */
    virtual const PlanarArcGeometry &getGeometry() const override {return planarArcGeometry_m;}

    /** Set LossDataSink to sink.
     *
     *  @param sink The LossDataSink. Ring takes ownership of memory
     *         allocated to sink
     */
    void setLossDataSink(LossDataSink* sink);

    /** Get pointer to lossDataSink.
     *
     *  Ring still owns the memory to which lossDataSink points.
     */
    PartBunchBase<double, 3>* getLossDataSink() const;

    /** Set RefPartBunchBase<double, 3> to bunch.
     *
     *  @param sink The Bunch. Ring borrows memory allocated to bunch.
     *
     *  Note for compliance with style guide and compatibility with parent two
     *  pointer to RefPartBunchBase<double, 3> are stored; this keeps them aligned
     */
    void setRefPartBunch(PartBunchBase<double, 3>* bunch);

    /** Get pointer to RefPartBunchBase<double, 3> from the bunch.
     *
     *  Ring does not own this memory (so neither does caller).
     */
    PartBunchBase<double, 3>* getRefPartBunch() const;

    /** Set the harmonic number for RF (number of bunches in the ring) */
    void setHarmonicNumber(double cyclHarm) {cyclHarm_m = cyclHarm;}

    /** Get the harmonic number for RF (number of bunches in the ring) */
    double getHarmonicNumber() {return cyclHarm_m;}
    // note this is not a const method to follow parent

    /** Set the nominal RF frequency */
    void setRFFreq(double rfFreq) {rfFreq_m = rfFreq;}

    /** Get the nominal RF frequency */
    double getRFFreq() const {return rfFreq_m;}

    /** Set the initial beam radius */
    void setBeamRInit(double rInit) {beamRInit_m = rInit;}

    /** Get the initial beam radius */
    double getBeamRInit() const {return beamRInit_m;}

    /** Set the initial beam azimuthal angle */
    void setBeamPhiInit(double phiInit) {beamPhiInit_m = phiInit;}

    /** Get the initial beam azimuthal angle */
    double getBeamPhiInit() const {return beamPhiInit_m;}

    /** Set the initial beam radial momentum */
    void setBeamPRInit(double pRInit) {beamPRInit_m = pRInit;}

    /** Get the initial beam radial momentum */
    double getBeamPRInit() const {return beamPRInit_m;}

    /** Set the initial element's radius */
    void setLatticeRInit(double rInit) {latticeRInit_m = rInit;}

    /** Get the initial element's radius */
    double getLatticeRInit() const {return latticeRInit_m;}

    /** Set the initial element's azimuthal angle */
    void setLatticePhiInit(double phiInit) {latticePhiInit_m = phiInit;}

    /** Get the initial  element's azimuthal angle */
    double getLatticePhiInit() const {return latticePhiInit_m;}

    /** Get the initial element's start position in cartesian coordinates */
    Vector_t getNextPosition() const;

    /** Get the initial element's start normal in cartesian coordinates */
    Vector_t getNextNormal() const;

    /** Set the first element's horizontal angle
     *
     *  Set the angle in the ring plane with respect to the tangent vector
     */
    void setLatticeThetaInit(double thetaInit) {latticeThetaInit_m = thetaInit;}

    /** Get the first element's horizontal angle
     *
     *  Get the angle in the ring plane with respect to the tangent vector
     */
    double getLatticeThetaInit() const {return latticeThetaInit_m;}

    /** Set the rotational symmetry of the ring (number of cells) */
    void setSymmetry(double symmetry) {symmetry_m = symmetry;}

    /** Set the scaling factor for the fields */
    void setScale(double scale) {scale_m = scale;}

    /** Get the rotational symmetry of the ring (number of cells) */
    double getSymmetry() const {return symmetry_m;}

    /** Set flag for closure checking */
    void setIsClosed(bool isClosed) {isClosed_m = isClosed;}

    /** Get flag for closure checking */
    double getIsClosed() const {return isClosed_m;}

    /** Set the ring aperture limits
     *  - minR; the minimum radius allowed during tracking. Throws if minR < 0
     *  - maxR; the maximum radius allowed during tracking. Throws if maxR < 0
     *  Sets willDoRingAperture_m to true. */
    void setRingAperture(double minR, double maxR);

    /** Get the ring minimum */
    double getRingMinR() const {return ::sqrt(minR2_m);}

    /** Get the ring maximum */
    double getRingMaxR() const {return ::sqrt(maxR2_m);}

    /** Lock the ring
     *
     *  Lock the ring; apply closure checks and symmetry properties as required.
     *  Impose rule that start must be before end and switch objects around if
     *  this is not the case. Sort by startPosition azimuthal angle.
     *
     *  Sets isLocked_m to true. New elements can no longer be added (as it may
     *  break the symmetry/bound checking)
     */
    void lockRing();

    /** Get the last section placed or NULL if no sections were placed yet */
    RingSection* getLastSectionPlaced() const;

    /** Get the list of sections at position pos */
    std::vector<RingSection*> getSectionsAt(const Vector_t& pos);

    /** Convert from a Vector3D to a Vector_t */
    static inline Vector_t convert(const Vector3D& vec);

    /** Convert from a Vector_t to a Vector3D */
    static inline Vector3D convert(const Vector_t& vec);

  private:
    // Force end to have azimuthal angle > start unless crossing phi = pi/-pi
    void resetAzimuths();

    // check for closure; throw an exception is ring is not closed within
    // tolerance; enforce closure to floating point precision
    void checkAndClose();

    // build a map that maps section to sections that overlap it
    void buildRingSections();

    void rotateToCyclCoordinates(Euclid3D& euclid3d) const;

    // predicate for sorting
    static bool sectionCompare(RingSection const* const sec1,
                               RingSection const* const sec2);

    /** Disabled */
    Ring();

    /** Disabled */
    Ring& operator=(const Ring& ring);

    void checkMidplane(Euclid3D delta) const;
    Rotation3D getRotationStartToEnd(Euclid3D delta) const;

    PlanarArcGeometry planarArcGeometry_m;

    // points to same location as RefPartBunch_m on the child bunch, but we
    // rename to keep in line with style guide
    //
    // Ring borrows this memory
    PartBunchBase<double, 3> *refPartBunch_m;

    // store for particles out of the aperture
    //
    // Ring owns this memory
    LossDataSink *lossDS_m;

    // initial position of the beam
    double beamRInit_m;
    double beamPRInit_m;
    double beamPhiInit_m;

    // position, orientation of the first lattice element
    double latticeRInit_m;
    double latticePhiInit_m;
    double latticeThetaInit_m;

    // aperture cut on the ring (before any field maps tracking)
    // note we store r^2
    bool willDoAperture_m = false;
    double minR2_m;
    double maxR2_m;

    // Ring is locked - new elements cannot be added
    bool isLocked_m;

    // Set to false to enable a non-circular ring (e.g. some weird spiral
    // geometry)
    bool isClosed_m;

    // number of cells/rotational symmetry of the ring
    int symmetry_m;

    double scale_m = 1.;

    // rf harmonic number
    double cyclHarm_m;

    // nominal rf frequency
    double rfFreq_m;

    // vector of RingSection sorted by phi (Component placement)
    double phiStep_m;
    std::vector<RingSectionList> ringSections_m;
    RingSectionList section_list_m;

    // tolerance on checking for geometry consistency
    static const double lengthTolerance_m;
    static const double angleTolerance_m;
};

Vector_t Ring::convert(const Vector3D& vec_3d) {
    return Vector_t(vec_3d(0), vec_3d(1), vec_3d(2));
}

Vector3D Ring::convert(const Vector_t& vec_t) {
    return Vector3D(vec_t[0], vec_t[1], vec_t[2]);
}

#endif //#ifndef RING_H
