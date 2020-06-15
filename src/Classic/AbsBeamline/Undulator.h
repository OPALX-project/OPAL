#ifndef CLASSIC_Undulator_HH
#define CLASSIC_Undulator_HH

// ------------------------------------------------------------------------
// $RCSfile: Undulator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Undulator
//   Defines the abstract interface for a drift space.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"

#ifdef OPAL_FEL
#include "mithra/fdtdSC.h"
#include "mithra/classes.h"
#endif


// Class Undulator
// ------------------------------------------------------------------------
/// Interface for drift space.
//  Class Undulator defines the abstract interface for a drift space.

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
    virtual void solve(MITHRA::FdTdSC& solver, MITHRA::Mesh& mesh, MITHRA::Bunch& bunch, MITHRA::Seed& seed);
#endif
  
    virtual void finalise();

    virtual bool bends() const;

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    //set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices); // Philippe was here

    //set number of slices for map tracking
    std::size_t getNSlices() const; // Philippe was here

    virtual void setK(double k);
    virtual double getK() const;
    virtual void setLambda(double lambda);
    virtual double getLambda() const;
    virtual void setNumPeriods(int np);
    virtual int getNumPeriods() const;
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

    /// Mithra file
    std::string fname_m;

    /// Mesh size
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
