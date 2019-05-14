#ifndef CLASSIC_RBend3D_HH
#define CLASSIC_RBend3D_HH

// ------------------------------------------------------------------------
// $RCSfile: RBend3D.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend3D
//   Defines the abstract interface for a solenoid magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:32 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/BendBase.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/BMultipoleField.h"

template <class T, unsigned Dim>
class PartBunchBase;
class Fieldmap;
class MeshData;

// Class RBend3D
// ------------------------------------------------------------------------
/// Interface for solenoids.
//  Class RBend3D defines the abstract interface for solenoid magnets.


class RBend3D: public BendBase {

public:

    /// Constructor with given name.
    explicit RBend3D(const std::string &name);

    RBend3D();
    RBend3D(const RBend3D &);
    virtual ~RBend3D();

    /** Inheritable copy constructor */
    ElementBase* clone() const override;

    /** Return the cell geometry */
    BGeometryBase& getGeometry() override;

    /** Return the cell geometry */
    const BGeometryBase& getGeometry() const override;

    /** Return a dummy (0.) field value (what is this for?) */
    EMField &getField() override;

    /** Return a dummy (0.) field value (what is this for?) */
    const EMField &getField() const override;

    /// Apply visitor to RBend3D.
    virtual void accept(BeamlineVisitor &) const override;

    virtual void addKR(int i, double t, Vector_t &K) override;

    virtual void addKT(int i, double t, Vector_t &K) override;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) override;

    virtual bool apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) override;

    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) override;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double &kineticEnergy) override;

    virtual void goOffline() override;

    //  Assign the field filename.
    void setFieldMapFN(std::string fn);

    void setFast(bool fast);

    bool getFast() const;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

    virtual bool isInside(const Vector_t &r) const override;

    MeshData getSurfaceMesh() const;

    virtual double getExitAngle() const;
private:
    double trackRefParticleThrough(double dt, bool print = false);

    //  std::string name;                 /**< The name of the object*/
    Fieldmap *myFieldmap_m;
    double fieldAmplitudeError_m;         /**< scale multiplier error*/

    double startField_m;                  /**< startingpoint of field, m*/
    double lengthField_m;

    bool fast_m;

    StraightGeometry geometry_m;

    BMultipoleField dummyField_m;

    // Not implemented.
    void operator=(const RBend3D &);
};

inline
double RBend3D::getExitAngle() const {
    return getBendAngle() - getEntranceAngle();
}

#endif // CLASSIC_RBend3D_HH