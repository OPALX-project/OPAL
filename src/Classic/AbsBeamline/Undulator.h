//
// Class Undulator
//   Defines all the methods used by the Undulator element.
//   The Undulator element uses a full wave solver from the
//   MITHRA library, see <https://github.com/aryafallahi/mithra/>.
//
// Copyright (c) 2020, Arnau Albà, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// Implemented as part of the MSc thesis
// "Start-to-End Modelling of the AWA Micro-Bunched Electron Cooling POP-Experiment"
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_Undulator_HH
#define CLASSIC_Undulator_HH

#include <string>
#include <vector>

#include "AbsBeamline/Component.h"

class Undulator : public Component {
public:
    /// Constructor with given name.
    explicit Undulator(const std::string& name);

    Undulator();
    Undulator(const Undulator& right);
    virtual ~Undulator();

    /// Apply visitor to Undulator.
    virtual void accept(BeamlineVisitor&) const;

    virtual void initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField);

    using Component::apply;
    void apply(PartBunchBase<double, 3>* itsBunch, CoordinateSystemTrafo const& refToLocalCSTrafo);

    virtual void finalise();

    virtual bool bends() const;

    virtual ElementType getType() const;

    virtual void getDimensions(double& zBegin, double& zEnd) const;

    void setK(double k);
    double getK() const;
    void setLambda(double lambda);
    double getLambda() const;
    void setNumPeriods(unsigned int np);
    unsigned int getNumPeriods() const;
    void setAngle(double theta);
    double getAngle() const;
    void setFilename(const std::string& fname);
    const std::string& getFilename() const;
    void setMeshLength(const std::vector<double>& ml);
    std::vector<double> getMeshLength() const;
    void setMeshResolution(const std::vector<double>& mr);
    std::vector<double> getMeshResolution() const;
    void setTruncationOrder(unsigned int trunOrder);
    unsigned int getTruncationOrder() const;
    void setTotalTime(double tt);
    double getTotalTime() const;
    void setDtBunch(double dtb);
    double getDtBunch() const;
    void setHasBeenSimulated(bool hbs);
    bool getHasBeenSimulated() const;

private:
    /// The undulator parameter
    double k_m;

    /// Undulator period
    double lambda_m;

    /// Number of periods
    unsigned int numPeriods_m;

    /// Polarisation angle of the undulator field
    double angle_m;

    /// Mithra file with output information
    std::string fname_m;

    /// Size of computational domain
    std::vector<double> meshLength_m;

    /// Mesh dx, dy, dz
    std::vector<double> meshResolution_m;

    /// First or second order absorbing boundary conditions
    unsigned int truncationOrder_m;

    /// Total time to run undulator
    double totalTime_m;

    /// Time step for the bunch position update
    double dtBunch_m;

    /// Boolean to indicate whether this undulator has already been simulated
    bool hasBeenSimulated_m;

    // Not implemented.
    void operator=(const Undulator&);
};

#endif  // CLASSIC_Undulator_HH
