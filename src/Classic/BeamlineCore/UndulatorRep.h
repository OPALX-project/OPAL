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
};

#endif // CLASSIC_UndulatorRep_HH
