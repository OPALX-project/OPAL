#ifndef CLASSIC_Degrader_HH
#define CLASSIC_Degrader_HH

// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Degrader
//   Defines the abstract interface for a beam Degrader.
//   *** MISSING *** Degrader interface is still incomplete.
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
#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include <vector>

class LossDataSink;

// Class Degrader
// ------------------------------------------------------------------------
/// Abstract collimator.
//  Class Degrader defines the abstract interface for a collimator.

class Degrader: public Component {

public:

    /// Plane selection.
    enum Plane {
        /// Monitor is off (inactive).
        OFF,
        /// Monitor acts on x-plane.
        X,
        /// Monitor acts on y-plane.
        Y,
        /// Monitor acts on both planes.
        XY
    };

    /// Constructor with given name.
    explicit Degrader(const std::string &name);

    Degrader();
    Degrader(const Degrader &rhs);
    virtual ~Degrader();

    /// Apply visitor to Degrader.
    virtual void accept(BeamlineVisitor &) const;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B);

    virtual bool applyToReferenceParticle(const Vector_t &R,
                                          const Vector_t &P,
                                          const double &t,
                                          Vector_t &E,
                                          Vector_t &B);

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField);

    virtual void initialise(PartBunchBase<double, 3> *bunch);

    virtual void finalise();

    virtual bool bends() const;

    virtual void goOnline(const double &kineticEnergy);

    virtual void goOffline();

    virtual ElementBase::ElementType getType() const;

    virtual void getDimensions(double &zBegin, double &zEnd) const;

    std::string  getDegraderShape(); // AAA

    void setOutputFN(std::string fn);
    std::string getOutputFN();

    virtual bool isInMaterial(double z);

private:

    // Not implemented.
    void operator=(const Degrader &);

    std::string filename_m;               /**< The name of the outputfile*/

    std::vector<double> PosX_m;
    std::vector<double> PosY_m;
    std::vector<double> PosZ_m;
    std::vector<double> MomentumX_m;
    std::vector<double> MomentumY_m;
    std::vector<double> MomentumZ_m;
    std::vector<double> time_m;
    std::vector<int> id_m;

    std::unique_ptr<LossDataSink> lossDs_m;
};

#endif // CLASSIC_Degrader_HH