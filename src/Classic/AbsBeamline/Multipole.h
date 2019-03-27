#ifndef CLASSIC_Multipole_HH
#define CLASSIC_Multipole_HH

// ------------------------------------------------------------------------
// $RCSfile: Multipole.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
// Description:
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/BMultipoleField.h"

template <class T, unsigned Dim>
class PartBunchBase;
class Fieldmap;

// Class Multipole
// ------------------------------------------------------------------------
/// Interface for general multipole.
//  Class Multipole defines the abstract interface for magnetic multipoles.
//  The order n of multipole components runs from 1 to N and is dynamically
//  adjusted. It is connected with the number of poles by the table
//
//  [tab 2 b]
//  [ROW]1[&]dipole[/ROW]
//  [ROW]2[&]quadrupole[/ROW]
//  [ROW]3[&]sextupole[/ROW]
//  [ROW]4[&]octupole[/ROW]
//  [ROW]5[&]decapole[/ROW]
//  [ROW]n[&]multipole with 2*n poles[/ROW]
//  [/TAB]
//  Units for multipole strengths are Teslas / m**(n-1).

class Multipole: public Component {

public:

    /// Constructor with given name.
    explicit Multipole(const std::string &name);

    Multipole();
    Multipole(const Multipole &);
    virtual ~Multipole();

    /// Apply visitor to Multipole.
    virtual void accept(BeamlineVisitor &) const override;


    /// Get multipole field.
    virtual BMultipoleField &getField() override = 0;

    /// Get multipole field. Version for const object.
    virtual const BMultipoleField &getField() const override = 0;

    /// Get normal component.
    //  Return the normal component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the return value is zero.
    double getNormalComponent(int n) const;

    /// Get skew component.
    //  Return the skew component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the return value is zero.
    double getSkewComponent(int n) const;

    /// Set normal component.
    //  Set the normal component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the component is created.
    void setNormalComponent(int, double);

    /// Set normal component.
    //  Set the normal component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the component is created.
    void setNormalComponent(int, double, double);

    /// Set skew component.
    //  Set the skew component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the component is created.
    void setSkewComponent(int, double);

    /// Set skew component.
    //  Set the skew component of order [b]n[/b] in T/m**(n-1).
    //  If [b]n[/b] is larger than the maximum order, the component is created.
    void setSkewComponent(int, double, double);

    size_t getMaxNormalComponentIndex() const;
    size_t getMaxSkewComponentIndex() const;
    
    //set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices);
    
    //set number of slices for map tracking
    std::size_t getNSlices() const;

    bool isFocusing(unsigned int component) const;

    /// Get geometry.
    virtual StraightGeometry &getGeometry() override = 0;

    /// Get geometry.
    virtual const StraightGeometry &getGeometry() const override = 0;

    virtual void addKR(int i, double t, Vector_t &K) override;

    virtual void addKT(int i, double t, Vector_t &K) override;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) override;

    virtual bool apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) override;

    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) override;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

    virtual bool isInside(const Vector_t &r) const override;
private:
    void computeField(Vector_t R, Vector_t &E, Vector_t &B);

    // Not implemented.
    void operator=(const Multipole &);
    std::vector<double> NormalComponents;
    std::vector<double> NormalComponentErrors;
    std::vector<double> SkewComponents;
    std::vector<double> SkewComponentErrors;
    int max_SkewComponent_m;
    int max_NormalComponent_m;
    std::size_t nSlices_m;
};

inline
void Multipole::setNormalComponent(int n, double v) {
    setNormalComponent(n, v, 0.0);
}

inline
void Multipole::setSkewComponent(int n, double v) {
    setSkewComponent(n, v, 0.0);
}

inline
size_t Multipole::getMaxNormalComponentIndex() const {
    return NormalComponents.size();
}

inline
size_t Multipole::getMaxSkewComponentIndex() const {
    return SkewComponents.size();
}

#endif // CLASSIC_Multipole_HH
