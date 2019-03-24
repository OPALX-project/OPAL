#ifndef CLASSIC_Solenoid_HH
#define CLASSIC_Solenoid_HH

// ------------------------------------------------------------------------
// $RCSfile: Solenoid.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Solenoid
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

#include "AbsBeamline/Component.h"

template <class T, unsigned Dim>
class PartBunchBase;
class Fieldmap;

// Class Solenoid
// ------------------------------------------------------------------------
/// Interface for solenoids.
//  Class Solenoid defines the abstract interface for solenoid magnets.


class Solenoid: public Component {

public:

    /// Constructor with given name.
    explicit Solenoid(const std::string &name);

    Solenoid();
    Solenoid(const Solenoid &);
    virtual ~Solenoid();

    /// Apply visitor to Solenoid.
    virtual void accept(BeamlineVisitor &) const;

    /// Get solenoid field Bz in Teslas.
    virtual double getBz() const = 0;

    void setKS(double ks);
    void setDKS(double ks);

    virtual void addKR(int i, double t, Vector_t &K);

    virtual void addKT(int i, double t, Vector_t &K);

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B);

    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B);

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    //  Assign the field filename.
    void setFieldMapFN(std::string fn);

    void setFast(bool fast);

    bool getFast() const;

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    virtual bool isInside(const Vector_t &r) const;

    virtual double getElementLength() const;

    virtual void getElementDimensions(double &zBegin, double &zEnd) const;

    virtual CoordinateSystemTrafo getEdgeToBegin() const;
    virtual CoordinateSystemTrafo getEdgeToEnd() const;

private:

    //  std::string name;                   /**< The name of the object*/
    std::string filename_m;               /**< The name of the inputfile*/
    Fieldmap *myFieldmap_m;
    double scale_m;                /**< scale multiplier*/
    double scaleError_m;                /**< scale multiplier error*/

    double startField_m;           /**< startingpoint of field, m*/
    double length_m;

    bool fast_m;
    // Not implemented.
    void operator=(const Solenoid &);
};



inline
CoordinateSystemTrafo Solenoid::getEdgeToBegin() const
{
    CoordinateSystemTrafo ret(Vector_t(0, 0, startField_m),
                              Quaternion(1, 0, 0, 0));

    return ret;
}

inline
CoordinateSystemTrafo Solenoid::getEdgeToEnd() const
{
    CoordinateSystemTrafo ret(Vector_t(0, 0, startField_m + length_m),
                              Quaternion(1, 0, 0, 0));

    return ret;
}
#endif // CLASSIC_Solenoid_HH