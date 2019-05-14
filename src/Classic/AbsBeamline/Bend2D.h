#ifndef CLASSIC_BEND_H
#define CLASSIC_BEND_H

// ------------------------------------------------------------------------
// $RCSfile: SBend.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: Bend2D
//   Defines the abstract interface for a general bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BendBase.h"
#include "Steppers/BorisPusher.h"
#include "Utilities/GeneralClassicException.h"

#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp.h"

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

#include <array>
#include <vector>

class Fieldmap;
class MeshData;

/*
 * Class Bend2D
 *
 * Interface for general bend magnet.
 *
 * ------------------------------------------------------------------------
 *
 */

class Bend2D: public BendBase {

public:

    /// Constructor with given name.
    explicit Bend2D(const std::string &name);

    Bend2D();
    Bend2D(const Bend2D &);
    virtual ~Bend2D();

    /// Apply visitor to Bend2D.
    virtual void accept(BeamlineVisitor &) const = 0;

    /*
     * Methods for OPAL-SLICE.
     */
    virtual void addKR(int i, double t, Vector_t &K) { };
    virtual void addKT(int i, double t, Vector_t &K) { };


    /*
     * Methods for OPAL-T.
     * ===================
     */

    /// Apply field to particles with coordinates in magnet frame.
    virtual bool apply(const size_t &i,
                       const double &t,
                       Vector_t &E,
                       Vector_t &B);

    /// Apply field to particles in beam frame.
    virtual bool apply(const Vector_t &R,
                       const Vector_t &P,
                       const double &t,
                       Vector_t &E,
                       Vector_t &B);

    virtual bool applyToReferenceParticle(const Vector_t &R,
                                          const Vector_t &P,
                                          const double &t,
                                          Vector_t &E,
                                          Vector_t &B);

    virtual void goOnline(const double &kineticEnergy);

    virtual void finalise();
    virtual void getDimensions(double &sBegin, double &sEnd) const;
    virtual ElementBase::ElementType getType() const = 0;
    virtual void initialise(PartBunchBase<double, 3> *bunch,
                            double &startField,
                            double &endField);

    double getBendRadius() const;
    double getEffectiveCenter() const;
    double getEffectiveLength() const;

    double getStartElement() const;

    /// Set quadrupole field component.
    void setK1(double k1);

    void resetReinitializeFlag();
    void resetRecalcRefTrajFlag();

    void setExitAngle(double exitAngle);
    virtual double getExitAngle() const;

    double getMapLength() const;

    std::vector<Vector_t> getOutline() const;
    MeshData getSurfaceMesh() const;

    virtual CoordinateSystemTrafo getEdgeToEnd() const;
    CoordinateSystemTrafo getBeginToEnd_local() const;

    virtual bool isInside(const Vector_t &r) const;


    //set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices);

    //set number of slices for map tracking
    std::size_t getNSlices() const;

    /// Get entrance fringe field length.
    //  Used to create fringe fields in ThickTracker, (before edge[m], after edge[m])
    std::array<double,2> getEntranceFringeFieldLength() const;

   /// Get exit fringe field length.
   //  Used to create fringe fields in ThickTracker, (before edge[m], after edge[m])
   std::array<double,2> getExitFringeFieldLength() const;

protected:
    void setMessageHeader(const std::string & header);
    double getStartField() const;
    Fieldmap* getFieldmap();

private:

#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(Maxwell, Zeros);
#endif
    // Not implemented.
    void operator=(const Bend2D &);

    void adjustFringeFields(double ratio);
    double calculateBendAngle();
    void calcEngeFunction(double zNormalized,
                          const std::vector<double> &engeCoeff,
                          int polyOrder,
                          double &engeFunc,
                          double &engeFuncDeriv,
                          double &engeFuncSecDerivNorm);
    Vector_t calcCentralField(const Vector_t &R,
                              double deltaX);
    Vector_t calcEntranceFringeField(const Vector_t &R,
                                     double deltaX);
    Vector_t calcExitFringeField(const Vector_t &R,
                                 double deltaX);
    void setupFringeWidths();
    bool calculateMapField(const Vector_t &R,
                           Vector_t &B);
    void calculateRefTrajectory(double &angleX,
                                double &angleY);
    double estimateFieldAdjustmentStep(double actualBendAngle,
                                       double mass,
                                       double betaGamma);
    void findBendEffectiveLength(double startField,
                                 double endField);
    void findBendStrength(double mass,
                          double gamma,
                          double betaGamma,
                          double charge);
    virtual bool findChordLength(Inform &msg,
                                 double &chordLength) = 0;
    bool findIdealBendParameters(double chordLength);
    void findReferenceExitOrigin(double &x, double &z);
    bool initializeFieldMap(Inform &msg);
    bool inMagnetCentralRegion(const Vector_t &R) const;
    bool inMagnetEntranceRegion(const Vector_t &R) const;
    bool inMagnetExitRegion(const Vector_t &R) const;
    bool isPositionInEntranceField(const Vector_t &R) const;
    bool isPositionInExitField(const Vector_t &R) const;
    void print(Inform &msg, double bendAngleX, double bendAngle);
    void readFieldMap(Inform &msg);
    bool reinitialize();
    void setBendEffectiveLength(double startField, double endField);
    void setBendStrength();
    void setEngeOriginDelta(double delta);
    void setFieldCalcParam();
    void setGapFromFieldMap();
    bool setupBendGeometry(Inform &msg, double &startField, double &endField);
    bool setupDefaultFieldMap(Inform &msg);
    void setFieldBoundaries(double startField, double endField);
    void setupPusher(PartBunchBase<double, 3> *bunch);
    bool treatAsDrift(Inform &msg, double chordlength);

    void setCSTrafoToEntranceRegion(const CoordinateSystemTrafo &trafo);
    void setCSTrafoToExitRegion(const CoordinateSystemTrafo &trafo);
    Vector_t transformToEntranceRegion(const Vector_t &R) const;
    Vector_t transformToExitRegion(const Vector_t &R) const;

    std::string messageHeader_m;

    BorisPusher pusher_m;       /// Pusher used to integrate reference particle
    /// through the bend.

    Fieldmap *fieldmap_m;       /// Magnet field map.
    bool fast_m;                /// Flag to turn on fast field calculation.

    double designRadius_m;      /// Bend design radius (m).

    /// and the entrance face of the magnet (radians).
    double exitAngle_m;         /// Angle between outgoing reference trajectory
    /// and the exit face of the magnet (radians).
    double fieldIndex_m;        /// Dipole field index.
    double startField_m;        /// Start of magnet field map in s coordinates (m).
    double endField_m;          /// End of magnet field map in s coordinates (m).

    double widthEntranceFringe_m;
    double widthExitFringe_m;

    /*
     * Flag to reinitialize the bend the first time the magnet
     * is called. This redefines the design energy of the bend
     * to the current average beam energy, keeping the bend angle
     * constant.
     */
    bool reinitialize_m;

    bool recalcRefTraj_m;       /// Re-calculate reference trajectory.

    /*
     * Enge function field map members.
     */

    /*
     * Entrance and exit position parameters. Ultimately they are used to
     * determine the origins of the entrance and exit edge Enge functions and
     * the extent of the field map. However, how they are used to do this
     * depends on how the bend using the map is setup in the OPAL input file.
     * So, we use generic terms to start.
     */
    double entranceParameter1_m;
    double entranceParameter2_m;
    double entranceParameter3_m;
    double exitParameter1_m;
    double exitParameter2_m;
    double exitParameter3_m;

    /// Enge coefficients for map entry and exit regions.
    std::vector<double> engeCoeffsEntry_m;
    std::vector<double> engeCoeffsExit_m;

    gsl_spline** entryFieldValues_m;
    gsl_spline** exitFieldValues_m;
    gsl_interp_accel *entryFieldAccel_m;
    gsl_interp_accel *exitFieldAccel_m;

    /*
     * All coordinates are with respect to (x, z) = (0, 0). It is
     * assumed the ideal reference trajectory passes through this point.
     */
    double deltaBeginEntry_m;       /// Perpendicular distance from entrance Enge
    /// function origin where Enge function starts.
    double deltaEndEntry_m;         /// Perpendicular distance from entrance Enge
    /// function origin that Enge function ends.
    int polyOrderEntry_m;           /// Enge function order for entry region.

    /*
     * The ideal reference trajectory passes through (xExit_m, zExit_m).
     */
    double xExit_m;
    double zExit_m;

    double deltaBeginExit_m;        /// Perpendicular distance from exit Enge
    /// function origin that Enge function starts.
    double deltaEndExit_m;          /// Perpendicular distance from exit Enge
    /// function origin that Enge function ends.
    int polyOrderExit_m;            /// Enge function order for entry region.

    double cosEntranceAngle_m;
    double sinEntranceAngle_m;
    double tanEntranceAngle_m;
    double tanExitAngle_m;

    CoordinateSystemTrafo beginToEnd_m;
    CoordinateSystemTrafo beginToEnd_lcs_m; // local coordinate system
    CoordinateSystemTrafo toEntranceRegion_m;
    CoordinateSystemTrafo toExitRegion_m;

    CoordinateSystemTrafo computeAngleTrafo_m;
    double maxAngle_m;

    std::size_t nSlices_m;
};


inline
void Bend2D::finalise() {
    online_m = false;
}

inline
void Bend2D::getDimensions(double &sBegin, double &sEnd) const {
    sBegin = startField_m;
    sEnd = endField_m;
}

inline
double Bend2D::getBendRadius() const {
    return designRadius_m;
}

inline
double Bend2D::getEffectiveCenter() const {
    return elementEdge_m + designRadius_m * angle_m / 2.0;
}

inline
double Bend2D::getEffectiveLength() const {
    return designRadius_m * angle_m;
}

inline
double Bend2D::getStartElement() const {
    return elementEdge_m;
}

inline
void Bend2D::setK1(double k1) {
    if (std::abs(k1) > 0.0) {
        throw GeneralClassicException("Bend2D::setK1",
                                      "Quadrupole field temporarily not supported");
    }
    fieldIndex_m = k1;
}

inline
void Bend2D::resetReinitializeFlag() {
    reinitialize_m = true;
}

inline
void Bend2D::resetRecalcRefTrajFlag() {
    recalcRefTraj_m = true;
}

inline
void Bend2D::setMessageHeader(const std::string & header)
{
    messageHeader_m = header;
}

inline
double Bend2D::getStartField() const
{
    return startField_m;
}

inline
Fieldmap* Bend2D::getFieldmap()
{
    return fieldmap_m;
}

inline
double Bend2D::getExitAngle() const
{
    return exitAngle_m;
}

inline
void Bend2D::setExitAngle(double angle)
{
    exitAngle_m = angle;
}

inline
double Bend2D::getMapLength() const
{
    return exitParameter2_m - entranceParameter2_m;
}

inline
CoordinateSystemTrafo Bend2D::getEdgeToEnd() const
{
    return beginToEnd_m;
}

inline
CoordinateSystemTrafo Bend2D::getBeginToEnd_local() const
{
    return beginToEnd_lcs_m;
}

inline
void Bend2D::setCSTrafoToEntranceRegion(const CoordinateSystemTrafo &trafo) {
    toEntranceRegion_m = trafo;
}

inline
void Bend2D::setCSTrafoToExitRegion(const CoordinateSystemTrafo &trafo) {
    toExitRegion_m = trafo;
}

inline
Vector_t Bend2D::transformToEntranceRegion(const Vector_t &R) const {
    return toEntranceRegion_m.transformTo(R);
}

inline
Vector_t Bend2D::transformToExitRegion(const Vector_t &R) const {
    return toExitRegion_m.transformTo(R);
}

#endif // CLASSIC_BEND_H