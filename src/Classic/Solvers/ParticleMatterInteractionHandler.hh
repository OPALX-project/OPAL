#ifndef PARTICLEMATTERINTERACTIONHANDLER_HH
#define PARTICLEMATTERINTERACTIONHANDLER_HH

#include <string>
#include "Algorithms/Vektor.h"

class ElementBase;

template <class T, unsigned Dim>
class PartBunchBase;
class Inform;

class ParticleMatterInteractionHandler {
public:
    ParticleMatterInteractionHandler(std::string name, ElementBase *elref);
    virtual ~ParticleMatterInteractionHandler() { };
    virtual void apply(PartBunchBase<double, 3> *bunch,
                       const std::pair<Vector_t, double> &boundingSphere,
                       size_t numParticlesInSimulation = 0) = 0;
    virtual const std::string getType() const = 0;
    virtual void print(Inform& os) = 0;
    virtual bool stillActive() = 0;
    virtual bool stillAlive(PartBunchBase<double, 3> *bunch) = 0;
    virtual double getTime() = 0;
    virtual std::string getName() = 0;
    virtual size_t getParticlesInMat() = 0;
    virtual unsigned getRediffused() = 0;
    virtual unsigned int getNumEntered() = 0;
    void setFlagAllParticlesIn(bool p);
    bool getFlagAllParticlesIn() const;
    void updateElement(ElementBase *newref);
    ElementBase* getElement();

    virtual bool computeEnergyLoss(Vector_t &P,//double &Eng,
                                   const double deltat,
                                   bool includeFluctuations = true) const = 0;
protected:
    ElementBase *element_ref_m;
    bool allParticleInMat_m; ///< if all particles are in matter stay inside the particle matter interaction
private:
    const std::string name_m;

};

inline
ParticleMatterInteractionHandler::ParticleMatterInteractionHandler(std::string name, ElementBase *elref):
    element_ref_m(elref),
    allParticleInMat_m(false),
    name_m(name)
{}

inline
void ParticleMatterInteractionHandler::updateElement(ElementBase *newref) {
    element_ref_m = newref;
}

inline
ElementBase* ParticleMatterInteractionHandler::getElement() {
    return element_ref_m;
}

inline
void ParticleMatterInteractionHandler::setFlagAllParticlesIn(bool p) {
  allParticleInMat_m = p;
}

inline
bool ParticleMatterInteractionHandler::getFlagAllParticlesIn() const {
    return allParticleInMat_m;
}
#endif // PARTICLEMATTERINTERACTION_HH