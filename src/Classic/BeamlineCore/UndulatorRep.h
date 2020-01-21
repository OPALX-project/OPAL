#ifndef CLASSIC_UndulatorRep_HH
#define CLASSIC_UndulatorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: UndulatorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: UndulatorRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Undulator.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/NullField.h"


// Class UndulatorRep
// ------------------------------------------------------------------------
/// Representation for a drift space.

class UndulatorRep: public Undulator {

public:

    /// Constructor with given name.
    explicit UndulatorRep(const std::string &name);

    UndulatorRep();
    UndulatorRep(const UndulatorRep &);
    virtual ~UndulatorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Get field.
    //  Version for non-constant object.
    virtual NullField &getField();

    /// Get field.
    //  Version for constant object.
    virtual const NullField &getField() const;

    /// Get geometry.
    //  Version for non-constant object.
    virtual StraightGeometry &getGeometry();

    /// Get geometry.
    //  Version for constant object.
    virtual const StraightGeometry &getGeometry() const;

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    virtual void setK(double k);
    virtual double getK() const;
    virtual void setLambda(double lambda);
    virtual double getLambda() const;
    virtual void setFilename(const std::string& fname);
    virtual const std::string& getFilename() const;
    virtual void setBF(double bf);
    virtual double getBF() const;
    virtual void setRadiationZ(std::vector<double> rz);
    virtual std::vector<double> getRadiationZ() const;
    virtual void setRadiationLambda(std::vector<double> rl);
    virtual std::vector<double> getRadiationLambda() const;
    virtual void setRadiationDirectory(const std::string& rd);
    virtual const std::string& getRadiationDirectory() const;
    virtual void setMeshLength(std::vector<double> ml);
    virtual std::vector<double> getMeshLength() const;
    virtual void setMeshResolution(std::vector<double> mr);
    virtual std::vector<double> getMeshResolution() const;
    virtual void setTruncationOrder(unsigned int trunOrder);
    virtual unsigned int getTruncationOrder() const;
    virtual void setSpaceCharge(bool sc);
    virtual bool getSpaceCharge() const;
    virtual void setEmitParticles(bool ep);
    virtual bool getEmitParticles() const;
    virtual void setTimeStepRatio(double m);
    virtual unsigned int getTimeStepRatio() const;
    virtual void setTotalTime(double tt);
    virtual double getTotalTime() const;
    virtual void setLFringe(double lf);
    virtual double getLFringe() const;
    virtual void setIsDone();
    virtual bool getIsDone() const;


private:

    // Not implemented.
    void operator=(const UndulatorRep &);

    /// The zero magnetic field.
    NullField field;

    /// The geometry.
    StraightGeometry geometry;

    /// The undulator parameter
    double k_m;

    /// Undulator period
    double lambda_m;

    /// Mithra file
    std::string fname_m;
    
    /// Desired bunching factor before entering undulator
    double bf_m;
    
    /// Distances from the bunch at which to measure radiation
    std::vector<double> radiationZ_m;

    /// The wavelength of the harmonic whose power should be measured, normalized to the radiation wavelength
    std::vector<double> radiationLambda_m;

    /// Directory in which to save radiation power
    std::string radiationDirectory_m;

    /// Mesh size
    std::vector<double> meshLength_m;

    /// Mesh dx, dy, dz
    std::vector<double> meshResolution_m;    

    /// First or second order absorbing boundary conditions
    unsigned int truncationOrder_m;

    /// Whether or not to consider space-charge effects
    bool spaceCharge_m;

    /// Emit particles or inject them all at once
    bool emitParticles_m;

    /// Ratio between mesh time-step and bunch time-step (both computed by mithra)
    unsigned int m_m;
    
    /// Total time to run undulator
    double totalTime_m;

    /// Distance from the undulator at which the bunch fields its fields
    double lFringe_m;
    
    /// Mithra has already run
    bool is_done_m;
};

#endif // CLASSIC_UndulatorRep_HH
