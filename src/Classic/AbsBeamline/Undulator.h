//
// Class Undulator
// Defines all the methods used by the Undulator element.
// The Undulator element uses a full wave solver from the
// the MITHRA library, see <https://github.com/aryafallahi/mithra/>.
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

#include "AbsBeamline/Component.h"

class Undulator: public Component {

public:
    /// Constructor with given name.
    explicit Undulator(const std::string &name);

    Undulator();
    Undulator(const Undulator &right);
    virtual ~Undulator();

    /// Apply visitor to Undulator.
    virtual void accept(BeamlineVisitor &) const;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

#ifdef OPAL_FEL
    virtual void apply(PartBunchBase<double, 3> *itsBunch, CoordinateSystemTrafo const& refToLocalCSTrafo);
#endif
  
    virtual void finalise();

    virtual bool bends() const;

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    //set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices);

    //set number of slices for map tracking
    std::size_t getNSlices() const;

    virtual void setK(double k);
    virtual double getK() const;
    virtual void setLambda(double lambda);
    virtual double getLambda() const;
    virtual void setNumPeriods(int np);
    virtual int getNumPeriods() const;
    virtual void setAngle(double theta);
    virtual double getAngle() const;
    virtual void setFilename(const std::string& fname);
    virtual const std::string& getFilename() const;
    virtual void setMeshLength(std::vector<double> ml);
    virtual std::vector<double> getMeshLength() const;
    virtual void setMeshResolution(std::vector<double> mr);
    virtual std::vector<double> getMeshResolution() const;
    virtual void setTruncationOrder(unsigned int trunOrder);
    virtual unsigned int getTruncationOrder() const;
    virtual void setTotalTime(double tt);
    virtual double getTotalTime() const;
    virtual void setDtBunch(double dtb);
    virtual double getDtBunch() const;
    virtual void setHasBeenSimlated(bool hbs);
    virtual bool getHasBeenSimulated() const;

private:

    double startField_m;
    std::size_t nSlices_m;
  
    /// The undulator parameter
    double k_m;

    /// Undulator period
    double lambda_m;

    /// Number of periods
    double numPeriods_m;

    /// Polarisation angle of the undulator field
    double angle_m;

    /// Length of the undulator's fringe field.
    double lFringe_m;

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
    void operator=(const Undulator &);
};

#endif // CLASSIC_Undulator_HH
