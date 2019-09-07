#ifndef CLASSIC_SOURCE_HH
#define CLASSIC_SOURCE_HH

#include "AbsBeamline/Component.h"
#include "Structure/LossDataSink.h"

class OpalBeamline;

template <class T, unsigned Dim>
class PartBunchBase;

class Source: public Component {

public:

    /// Constructor with given name.
    explicit Source(const std::string &name);

    Source();
    Source(const Source &);
    virtual ~Source();

    /// Apply visitor to Source.
    virtual void accept(BeamlineVisitor &) const override;

    virtual void addKR(int i, double t, Vector_t &K) override;

    virtual void addKT(int i, double t, Vector_t &K) override;

    virtual bool apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) override;

    virtual void initialise(PartBunchBase<double, 3> *bunch, double &startField, double &endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double &kineticEnergy) override;

    virtual void goOffline() override;

    virtual ElementBase::ElementType getType() const override;

    virtual void getDimensions(double &zBegin, double &zEnd) const override;

private:

    double ElementEdge_m;
    double startField_m;           /**< startingpoint of field, m*/
    double endField_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    // Not implemented.
    void operator=(const Source &);
};
#endif // CLASSIC_SOURCE_HH