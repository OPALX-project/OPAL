#ifndef CLASSIC_FlexibleCollimator_HH
#define CLASSIC_FlexibleCollimator_HH

#include "AbsBeamline/Component.h"
#include "Utilities/MSLang.h"
#include "Utilities/MSLang/QuadTree.h"

#include <memory>

class BeamlineVisitor;
class LossDataSink;
// Class FlexibleCollimator
// ------------------------------------------------------------------------
/// Abstract collimator.
//  Class FlexibleCollimator defines the abstract interface for a collimator.

class FlexibleCollimator: public Component {

public:

    /// Constructor with given name.
    explicit FlexibleCollimator(const std::string &name);

    FlexibleCollimator();
    FlexibleCollimator(const FlexibleCollimator &rhs);
    virtual ~FlexibleCollimator();

    /// Apply visitor to FlexibleCollimator.
    virtual void accept(BeamlineVisitor &) const override;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) override;

    virtual bool applyToReferenceParticle(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) override;

    virtual bool checkCollimator(PartBunchBase<double, 3> *bunch, const int turnnumber, const double t, const double tstep);

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void initialise(PartBunchBase<double, 3> *bunch);

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double &kineticEnergy) override;

    virtual void goOffline() override;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

    void print();

    void setOutputFN(const std::string &fn);
    std::string getOutputFN() const;

    unsigned int getLosses() const;

    void setDescription(const std::string &desc);
    std::string getDescription() const;

    bool isStopped(const Vector_t &R, const Vector_t &P, double recpgamma);

    void writeHolesAndQuadtree(const std::string &baseFilename) const;

private:

    // Not implemented.
    void operator=(const FlexibleCollimator &);

    std::string description_m;
    std::vector<std::shared_ptr<mslang::Base> > holes_m;
    mslang::BoundingBox bb_m;
    mslang::QuadTree tree_m;

    std::string filename_m;               /**< The name of the outputfile*/

    bool informed_m;
    unsigned int losses_m;
    std::unique_ptr<LossDataSink> lossDs_m;

    ParticleMatterInteractionHandler *parmatint_m;
};

inline
unsigned int FlexibleCollimator::getLosses() const {
    return losses_m;
}

inline
bool FlexibleCollimator::bends() const {
    return false;
}

inline
void FlexibleCollimator::setOutputFN(const std::string &fn) {
    filename_m = fn;
}

inline
std::string FlexibleCollimator::getOutputFN() const {
    return filename_m;
}

inline
std::string FlexibleCollimator::getDescription() const {
    return description_m;
}

#endif // CLASSIC_FlexibleCollimator_HH
