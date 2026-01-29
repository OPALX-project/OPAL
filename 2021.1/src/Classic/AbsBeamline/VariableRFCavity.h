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

#ifndef CLASSIC_ABSBEAMLINE_VariableRFCavity_HH
#define CLASSIC_ABSBEAMLINE_VariableRFCavity_HH

#include "Algorithms/AbstractTimeDependence.h"
#include "Fields/EMField.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "AbsBeamline/Component.h"

class Fieldmap;

/** @class VariableRFCavity
 *
 *  Generates a field like
 *      E = E0*a(t)*sin{f(t)*t-q(t)}
 *      B = 0
 *  where E0, B0 are user defined fields, a(t), f(t), q(t) are time
 *  dependent amplitude, frequency, phase respectively; it is assumed that these
 *  quantities vary sufficiently slowly that Maxwell is satisfied.
 * 
 *  The time dependent quantities are 
 */
class VariableRFCavity: public Component {
  public:
    /// Constructor with given name.
    explicit VariableRFCavity(const std::string &name);
    /** Copy Constructor; performs deepcopy on time-dependence models */
    VariableRFCavity(const VariableRFCavity &);
    /** Default constructor */
    VariableRFCavity();
    /** Assignment operator; performs deepcopy on time-dependence models*/
    VariableRFCavity& operator=(const VariableRFCavity &);
    /** Destructor does nothing
     * 
     * The shared_ptrs will self-destruct when reference count goes to 0
     */
    virtual ~VariableRFCavity();

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

    /** Initialise ready for tracking
     * 
     *  Just sets RefPartBunch_m
     */
    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    /** Finalise following tracking
     * 
     *  Just sets RefPartBunch_m to NULL
     */
    virtual void finalise() override;

    /** @returns false (cavity does not bend the trajectory) */
    virtual bool bends() const override {return false;}

    /** Not used (does nothing) */
    virtual void getDimensions(double &/*zBegin*/, double &/*zEnd*/) const override {}

    /** Get the amplitude at a given time
     *
     *  @param time: the time at which the amplitude is calculated
     *
     *  @returns the RF field gradient.
     */
    virtual inline double getAmplitude(double time) const;

    /** Get the frequency at a given time
     *
     *  @param time: the time at which the frequency is calculated
     *
     *  @returns the RF cavity frequency.
     */
    virtual inline double getFrequency(double time) const;

    /** Get the phase at a given time
     *
     *  @param phase: the time at which the phase is calculated
     *
     *  @returns the RF cavity phase.
     */
    virtual inline double getPhase(double time) const;

    /** @returns the full height of the cavity */
    virtual inline double getHeight() const;
    /** @returns the full width of the cavity */
    virtual inline double getWidth() const;
    /** @returns the length of the cavity */
    virtual inline double getLength() const;
    /** Set the full height of the cavity */
    virtual inline void setHeight(double fullHeight);
    /** Set the full width of the cavity */
    virtual inline void setWidth(double fullWidth);
    /** Set the length of the cavity */
    virtual void setLength(double length);

    /** @returns shared_ptr to the amplitude (field gradient) time dependence */
    virtual std::shared_ptr<AbstractTimeDependence> getAmplitudeModel() const;
    /** @returns shared_ptr to the phase time dependence */
    virtual std::shared_ptr<AbstractTimeDependence> getPhaseModel() const;
    /** @returns shared_ptr to the frequency */
    virtual std::shared_ptr<AbstractTimeDependence> getFrequencyModel() const;

    /** Set the amplitude (field gradient) time dependence */
    virtual void setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> time_dep);
    /** Set the phase time dependence */
    virtual void setPhaseModel(std::shared_ptr<AbstractTimeDependence> time_dep);
    /** Set the frequency time dependence */
    virtual void setFrequencyModel(std::shared_ptr<AbstractTimeDependence> time_dep);

    /** Set the amplitude time dependence name
     *
     *  The name is used to find the amplitude model at accept time
     */
    virtual void setAmplitudeName(std::string amplitude)
    { amplitudeName_m = amplitude; }

    /** Set the phase time dependence name
     *
     *  The name is used to find the phase model at accept time
     */
    virtual void setPhaseName(std::string phase)
    { phaseName_m = phase; }

    /** Set the frequency time dependence name
     *
     *  The name is used to find the frequency model at accept time
     */
    virtual void setFrequencyName(std::string frequency)
    { frequencyName_m = frequency; }

    /** Set the cavity geometry */
    virtual StraightGeometry& getGeometry() override;
    /** @returns the cavity geometry */
    virtual const StraightGeometry& getGeometry() const override;

    /// Not implemented
    virtual EMField &getField() override;
    /// Not implemented
    virtual const EMField &getField() const override;
  protected:
    void initNull();
    void initialise() const;
    std::shared_ptr<AbstractTimeDependence> phaseTD_m;
    std::shared_ptr<AbstractTimeDependence> amplitudeTD_m;
    std::shared_ptr<AbstractTimeDependence> frequencyTD_m;
    std::string phaseName_m;
    std::string amplitudeName_m;
    std::string frequencyName_m;
    double halfWidth_m;
    double halfHeight_m;
    double _length;
    /// The cavity's geometry.
    StraightGeometry geometry;
    static const double lengthUnit_m;

private:
};

double VariableRFCavity::getAmplitude(double time) const {
    return amplitudeTD_m->getValue(time);
}

double VariableRFCavity::getPhase(double time) const {
    return phaseTD_m->getValue(time);
}

double VariableRFCavity::getFrequency(double time) const {
    return frequencyTD_m->getValue(time);
}

double VariableRFCavity::getHeight() const {
    return halfHeight_m*2/lengthUnit_m;
}

double VariableRFCavity::getWidth() const {
    return halfWidth_m*2/lengthUnit_m;
}

double VariableRFCavity::getLength() const {
    return _length/lengthUnit_m;
}
 
void VariableRFCavity::setHeight(double fullHeight) {
    halfHeight_m = fullHeight/2*lengthUnit_m;
}

void VariableRFCavity::setWidth(double fullWidth) {
    halfWidth_m = fullWidth/2*lengthUnit_m;
}




#endif // CLASSIC_VirtualRFCavity_HH
