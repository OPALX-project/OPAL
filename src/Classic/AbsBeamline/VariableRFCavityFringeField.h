/*
 *  Copyright (c) 2014, Chris Rogers
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

#ifndef CLASSIC_ABSBEAMLINE_VariableRFCavityFringeField_HH
#define CLASSIC_ABSBEAMLINE_VariableRFCavityFringeField_HH

#include <iostream>
#include "AbsBeamline/VariableRFCavity.h"

class Fieldmap;
namespace endfieldmodel {
    class EndFieldModel;
}

/** @class VariableRFCavityFringeField
 *
 *  Generates a field like
 *      Ey = E0*a(t)*y^{2n+1} g_n(z) sin{f(t)*t-q(t)}
 *      Ez = E0*a(t)*y^{2n}   f_n(z) sin{f(t)*t-q(t)}
 *      Bx = B0*a(t)*y^{2n+1} h_n(z) cos{f(t)*t-q(t)}
 *  where E0, is a user-defined field and B0 is the corresponding magnetic field
 *  magnitude; f_n is a user-defined axial field dependence and g_n, h_n are the
 *  corresponding off-axis field dependencies. a(t), f(t), q(t) are time
 *  dependent amplitude, frequency, phase respectively; it is assumed that these
 *  quantities vary sufficiently slowly that Maxwell is satisfied.
 *
 *  Inherits from VariableRFCavity; inheritance is used here to share code from
 *  VariableRFCavity.
 *
 *  Set/get methods use metres; but internally we store as mm (for tracking)
 * 
 *  Field units are kG and GV/mm
 */
class VariableRFCavityFringeField : public VariableRFCavity {
  public:
    /// Constructor with given name.
    explicit VariableRFCavityFringeField(const std::string &name);
    /** Copy Constructor; performs deepcopy on time-dependence models */
    VariableRFCavityFringeField(const VariableRFCavityFringeField &);
    /** Default constructor */
    VariableRFCavityFringeField();
    /** Assignment operator; performs deepcopy on time-dependence models*/
    VariableRFCavityFringeField& operator=(const VariableRFCavityFringeField &);
    /** Destructor does nothing
     * 
     * The shared_ptrs will self-destruct when reference count goes to 0
     */
    virtual ~VariableRFCavityFringeField();

    /** Apply visitor to RFCavity.
     *
     *  The RF cavity finds the "time dependence" models by doing a string
     *  lookup against a list held by AbstractTimeDependence at accept time.
     */
    virtual void accept(BeamlineVisitor &) const override;

    /** Inheritable deepcopy method */
    virtual ElementBase* clone() const override;

    /** Calculate the field at the position of the i^th particle 
     *
     *  @param i indexes the particle whose field we need
     *  @param t the time at which the field is calculated
     *  @param E return value with electric field strength
     *  @param B return value with magnetic field strength
     *
     *  @returns True if particle is outside the boundaries; else False
     */
    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) override;

    /** Calculate the field at a given position
     *
     *  @param R the position at which the field is calculated
     *  @param P the momentum (not used)
     *  @param t the time at which the field is calculated
     *  @param E return value; filled with electric field strength
     *  @param B return value; filled with magnetic field strength
     *
     *  @returns True if particle is outside the boundaries; else False
     */
    virtual bool apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) override;


    /** Calculate the field at a given position. This is identical to "apply".
     *
     *  @param R the position at which the field is calculated
     *  @param P the momentum (not used)
     *  @param t the time at which the field is calculated
     *  @param E return value; filled with electric field strength
     *  @param B return value; filled with magnetic field strength
     *
     *  @returns True if particle is outside the boundaries; else False
     */
    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) override;

    /** Set the endFieldModel
     *
     *  @param endField the new endFieldModel. VariableRFCavityFringe takes
     *  ownership of the memory allocated to endField.
     */
    virtual void setEndField
                       (std::shared_ptr<endfieldmodel::EndFieldModel> endField);

    /** Get the endFieldModel
     * 
     *  @returns the endFieldModel; VariableRFCavityFringe retains ownership of
     *  the memory allocated to the endFieldModel.
     */
    virtual inline std::shared_ptr<endfieldmodel::EndFieldModel>
                                                            getEndField() const;

    /** Initialise ready for tracking
     * 
     *  Does any setup on VirtualRFCavity and sets field expansion coefficients
     */
    virtual void initialise(PartBunchBase<double, 3> *bunch,
                            double &startField,
                            double &endField) override;

    /** Get the offset of centre of the cavity field from the element start in metres */
    virtual inline void setCavityCentre(double zCentre);
    /** Set the offset of centre of the cavity field from the element start in metres */
    virtual inline double getCavityCentre() const;
    /** Set the maximum y power */
    virtual inline void setMaxOrder(size_t maxOrder);
    /** Get the maximum y power that will be used in field calculations */
    virtual inline size_t getMaxOrder() const;
    /** Set the coefficients for calculating the field expansion */
    void initialiseCoefficients();
    /** Print the coefficients to ostream out */
    void printCoefficients(std::ostream& out) const;
    /** Get the coefficients for Ez */
    inline std::vector<std::vector<double> > getEzCoefficients() const;
    /** Get the coefficients for Ey */
    inline std::vector<std::vector<double> > getEyCoefficients() const;
    /** Get the coefficients for Bx */
    inline std::vector<std::vector<double> > getBxCoefficients() const;
protected:
    double zCentre_m; // offsets this element
    size_t maxOrder_m;
    std::shared_ptr<endfieldmodel::EndFieldModel> endField_m;
    std::vector<std::vector<double> > f_m;
    std::vector<std::vector<double> > g_m;
    std::vector<std::vector<double> > h_m;
private:
    void initNull();
};

void VariableRFCavityFringeField::setCavityCentre(double zCentre) {
    zCentre_m = zCentre*1000.; // stored internally as mm
}

double VariableRFCavityFringeField::getCavityCentre() const {
    return zCentre_m/1000.; // stored internally as mm
}

void VariableRFCavityFringeField::setMaxOrder(size_t maxOrder) {
    maxOrder_m = maxOrder;
    initialiseCoefficients();
}

size_t VariableRFCavityFringeField::getMaxOrder() const {
    return maxOrder_m;
}

std::shared_ptr<endfieldmodel::EndFieldModel>
                              VariableRFCavityFringeField::getEndField() const {
    return endField_m;
}

std::vector<std::vector<double> >
                        VariableRFCavityFringeField::getEzCoefficients() const {
    return f_m;
}

std::vector<std::vector<double> >
                        VariableRFCavityFringeField::getEyCoefficients() const {
    return g_m;
}

std::vector<std::vector<double> >
                        VariableRFCavityFringeField::getBxCoefficients() const {
    return h_m;
}

#endif // #ifdef CLASSIC_VirtualRFCavityFringeField_HH
