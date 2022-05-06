//
// Class Source
//   Defines the abstract interface for a source.
//
// Copyright (c) 200x - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_SOURCE_HH
#define CLASSIC_SOURCE_HH

#include "AbsBeamline/Component.h"
#include "Structure/LossDataSink.h"

class OpalBeamline;
class LossDataSink;

template <class T, unsigned Dim>
class PartBunchBase;

class Source: public Component {

public:

    /// Constructor with given name.
    explicit Source(const std::string& name);

    Source();
    Source(const Source&);
    virtual ~Source();

    /// Apply visitor to Source.
    virtual void accept(BeamlineVisitor&) const override;

    virtual bool apply(const size_t& i, const double& t, Vector_t& E, Vector_t& B) override;

    virtual void initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double& kineticEnergy) override;

    virtual void goOffline() override;

    virtual ElementType getType() const override;

    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    virtual int getRequiredNumberOfTimeSteps() const override;

    void setTransparent();


private:

    double startField_m; /**< startingpoint of field, m*/
    double endField_m;

    bool isTransparent_m;

    std::unique_ptr<LossDataSink> lossDs_m; /**< Handling for store the particle out of region*/

    // Not implemented.
    void operator=(const Source&);
};

inline
int Source::getRequiredNumberOfTimeSteps() const {
    return 0;
}
#endif // CLASSIC_SOURCE_HH
