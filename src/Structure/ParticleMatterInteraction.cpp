// ------------------------------------------------------------------------
// $RCSfile: ParticleMatterInteraction.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: ParticleMatterInteraction
//   The class for the OPAL PARTICLEMATTERINTERACTION command.
//
// $Date: 2009/07/14 22:09:00 $
// $Author: C. Kraus $
//
// ------------------------------------------------------------------------

#include "Structure/ParticleMatterInteraction.h"
#include "Solvers/CollimatorPhysics.hh"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"
#include "AbsBeamline/ElementBase.h"
#include "Utilities/Util.h"

extern Inform *gmsg;

using namespace Physics;


// Class ParticleMatterInteraction
// ------------------------------------------------------------------------

// The attributes of class ParticleMatterInteraction.
namespace {
    enum {
        // DESCRIPTION OF SINGLE PARTICLE:
        TYPE,       // The type of the wake
        MATERIAL,   // From of the tube
        ENABLERUTHERFORD,
        SIZE
    };
}

ParticleMatterInteraction::ParticleMatterInteraction():
    Definition(SIZE, "PARTICLEMATTERINTERACTION",
               "The \"SURFACE_PHYSICS\" statement defines data for the particle mater interaction handler "
               "on an element."),
    handler_m(0) {
    itsAttr[TYPE] = Attributes::makeString
                    ("TYPE", "Specifies the particle mater interaction handler: Collimator");

    itsAttr[MATERIAL] = Attributes::makeString
                        ("MATERIAL", "The material of the surface");

    itsAttr[ENABLERUTHERFORD] = Attributes::makeBool("ENABLERUTHERFORD",
                                                     "Enable large angle scattering", true);

    ParticleMatterInteraction *defParticleMatterInteraction = clone("UNNAMED_PARTICLEMATTERINTERACTION");
    defParticleMatterInteraction->builtin = true;

    try {
        defParticleMatterInteraction->update();
        OpalData::getInstance()->define(defParticleMatterInteraction);
    } catch(...) {
        delete defParticleMatterInteraction;
    }

    registerOwnership(AttributeHandler::STATEMENT);
}


ParticleMatterInteraction::ParticleMatterInteraction(const std::string &name, ParticleMatterInteraction *parent):
    Definition(name, parent),
    handler_m(parent->handler_m)
{}


ParticleMatterInteraction::~ParticleMatterInteraction() {
    if(handler_m)
        delete handler_m;
}


bool ParticleMatterInteraction::canReplaceBy(Object *object) {
    // Can replace only by another PARTICLEMATTERINTERACTION.
    return dynamic_cast<ParticleMatterInteraction *>(object) != 0;
}


ParticleMatterInteraction *ParticleMatterInteraction::clone(const std::string &name) {
    return new ParticleMatterInteraction(name, this);
}


void ParticleMatterInteraction::execute() {
    update();
}


ParticleMatterInteraction *ParticleMatterInteraction::find(const std::string &name) {
    ParticleMatterInteraction *parmatint = dynamic_cast<ParticleMatterInteraction *>(OpalData::getInstance()->find(name));

    if(parmatint == 0) {
        throw OpalException("ParticleMatterInteraction::find()", "ParticleMatterInteraction \"" + name + "\" not found.");
    }
    return parmatint;
}


void ParticleMatterInteraction::update() {
    // Set default name.
    if(getOpalName().empty()) setOpalName("UNNAMED_PARTICLEMATTERINTERACTION");
}


void ParticleMatterInteraction::initParticleMatterInteractionHandler(ElementBase &element) {
    *gmsg << "* ************* P A R T I C L E  M A T E R  I N T E R A C T I O N ****************** " << endl;
    *gmsg << "* ParticleMatterInteraction::initParticleMatterInteractionHandler " << endl;
    *gmsg << "* ********************************************************************************** " << endl;

    std::string material = Util::toUpper(Attributes::getString(itsAttr[MATERIAL]));
    bool enableRutherford = Attributes::getBool(itsAttr[ENABLERUTHERFORD]);

    const std::string type = Util::toUpper(Attributes::getString(itsAttr[TYPE]));
    if(type == "CCOLLIMATOR" ||
       type == "COLLIMATOR" ||
       type == "DEGRADER") {

        handler_m = new CollimatorPhysics(getOpalName(), &element, material, enableRutherford);
        *gmsg << *this << endl;
    } else {
        handler_m = 0;
        INFOMSG(getOpalName() + ": no particle mater interaction handler attached, TYPE == " << Attributes::getString(itsAttr[TYPE]) << endl);
    }

}

void ParticleMatterInteraction::updateElement(ElementBase *element) {
    handler_m->updateElement(element);
}

void ParticleMatterInteraction::print(std::ostream &os) const {
    os << "* ************* P A R T I C L E  M A T T E R  I N T E R A C T I O N ****************** " << std::endl;
    os << "* PARTICLEMATTERINTERACTION " << getOpalName() << '\n'
       << "* MATERIAL       " << Attributes::getString(itsAttr[MATERIAL]) << '\n';
    os << "* ********************************************************************************** " << std::endl;
}