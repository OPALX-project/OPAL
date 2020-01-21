// ------------------------------------------------------------------------
// $RCSfile: OpalUndulator.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalUndulator
//   The class of OPAL drift spaces.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalUndulator.h"
#include "Structure/BoundaryGeometry.h"
#include "Attributes/Attributes.h"
#include "BeamlineCore/UndulatorRep.h"
#include "Structure/OpalWake.h"
#include "Structure/ParticleMatterInteraction.h"

// Class OpalUndulator
// ------------------------------------------------------------------------

OpalUndulator::OpalUndulator():
    OpalElement(SIZE, "UNDULATOR",
                "The \"UNDULATOR\" element defines a undulator."),
    owk_m(NULL),
    parmatint_m(NULL),
    obgeo_m(NULL) {
    // CKR: the following 3 lines are redundant: OpalElement does this already!
    //      they prevent drift from working properly
    //
    //     itsAttr[LENGTH] = Attributes::makeReal
    //         ("LENGTH", "Undulator length");

    //     registerRealAttribute("LENGTH");
    itsAttr[GEOMETRY] = Attributes::makeString
                        ("GEOMETRY", "BoundaryGeometry for Undulators");

    itsAttr[NSLICES] = Attributes::makeReal
                          ("NSLICES",
                          "The number of slices/ steps for this element in Map Tracking", 1);

    itsAttr[K] = Attributes::makeReal
                          ("K",
                          "The undulator parameter", 1);

    itsAttr[LAMBDA] = Attributes::makeReal
                          ("LAMBDA",
                          "The undulator period", 0.0);

    itsAttr[FNAME] = Attributes::makeString
      ("FNAME", "jobfile with desired output data from the undulator", "");

    itsAttr[BF] = Attributes::makeReal
                          ("BF",
                          "Desired initial bunching factor. If zero, the bunch's density modulation is not modified", 0.0);
    itsAttr[RADZ] = Attributes::makeRealArray
                          ("RADZ",
			   "Distances from bunch in which to measure radiation");
    itsAttr[RADLAMBDA] = Attributes::makeRealArray
                          ("RADLAMBDA",
			   "Normalized wavelengths to measure");
    itsAttr[RADDIRECTORY] = Attributes::makeString
                          ("RADDIRECTORY",
                          "Directory to store radiation power", "");
    itsAttr[MESHLENGTH] = Attributes::makeRealArray
                          ("MESHLENGTH",
			   "Size of computational mesh");
    itsAttr[MESHRESOLUTION] = Attributes::makeRealArray
                          ("MESHRESOLUTION",
			   "dx, dy, dz of the mesh");
    itsAttr[TRUNORDER] = Attributes::makeReal
                          ("TRUNORDER",
                          "Order of boundary absorbing conditions", 2);
    itsAttr[SPACECHARGE] = Attributes::makeReal
                          ("SPACECHARGE",
                          "Whether to consider space-charge effects or not", 1);
    itsAttr[EMITPARTICLES] = Attributes::makeReal
                          ("EMITPARTICLES",
                          "Emit particles or inject the whole bunch at once", 1);
    itsAttr[TIMESTEPRATIO] = Attributes::makeReal
                          ("TIMESTEPRATIO",
                          "Ratio between mesh time-step and bunch time-step", 1);
    itsAttr[TOTALTIME] = Attributes::makeReal
                          ("TOTALTIME",
                          "Total time of undulator simulation", 0.0);
    itsAttr[LFRINGE] = Attributes::makeReal("LFRINGE", 
					    "Distance from the undulator where the bunch starts to feel the undulator fields", 1.0);


    registerStringAttribute("GEOMETRY");
    registerRealAttribute("NSLICES");
    registerRealAttribute("K");
    registerRealAttribute("LAMBDA");
    registerStringAttribute("FNAME");
    registerRealAttribute("BF");
    registerStringAttribute("RADDIRECTORY");
    registerRealAttribute("TRUNORDER");
    registerRealAttribute("SPACECHARGE");
    registerRealAttribute("EMITPARTICLES");
    registerRealAttribute("TIMESTEPRATIO");
    registerRealAttribute("TOTALTIME");
    registerRealAttribute("LFRINGE");
    registerOwnership();

    setElement(new UndulatorRep("UNDULATOR"));
}


OpalUndulator::OpalUndulator(const std::string &name, OpalUndulator *parent):
    OpalElement(name, parent),
    owk_m(NULL),
    parmatint_m(NULL),
    obgeo_m(NULL) {
    setElement(new UndulatorRep(name));
}


OpalUndulator::~OpalUndulator() {
    if(owk_m)
        delete owk_m;
    if(parmatint_m)
        delete parmatint_m;
    if(obgeo_m)
	delete obgeo_m;
}


OpalUndulator *OpalUndulator::clone(const std::string &name) {
    return new OpalUndulator(name, this);
}


bool OpalUndulator::isUndulator() const {
    return true;
}


void OpalUndulator::update() {
    OpalElement::update();

    UndulatorRep *ur = static_cast<UndulatorRep *>(getElement());
    ur->setElementLength(Attributes::getReal(itsAttr[LENGTH]));
    ur->setNSlices(Attributes::getReal(itsAttr[NSLICES]));
    if(itsAttr[WAKEF] && owk_m == NULL) {
        owk_m = (OpalWake::find(Attributes::getString(itsAttr[WAKEF])))->clone(getOpalName() + std::string("_wake"));
        owk_m->initWakefunction(*ur);
        ur->setWake(owk_m->wf_m);
    }

    if(itsAttr[PARTICLEMATTERINTERACTION] && parmatint_m == NULL) {
        parmatint_m = (ParticleMatterInteraction::find(Attributes::getString(itsAttr[PARTICLEMATTERINTERACTION])))->clone(getOpalName() + std::string("_parmatint"));
        parmatint_m->initParticleMatterInteractionHandler(*ur);
        ur->setParticleMatterInteraction(parmatint_m->handler_m);
    }
    if(itsAttr[GEOMETRY] && obgeo_m == NULL) {
        obgeo_m = (BoundaryGeometry::find(Attributes::getString(itsAttr[GEOMETRY])))->clone(getOpalName() + std::string("_geometry"));
        if(obgeo_m) {
	    ur->setBoundaryGeometry(obgeo_m);
        }
    }


    ur->setK(Attributes::getReal(itsAttr[K]));
    ur->setLambda(Attributes::getReal(itsAttr[LAMBDA]));
    ur->setFilename(Attributes::getString(itsAttr[FNAME]));
    ur->setBF(Attributes::getReal(itsAttr[BF]));
    ur->setRadiationZ(Attributes::getRealArray(itsAttr[RADZ]));
    ur->setRadiationLambda(Attributes::getRealArray(itsAttr[RADLAMBDA]));
    ur->setRadiationDirectory(Attributes::getString(itsAttr[RADDIRECTORY]));
    ur->setMeshLength(Attributes::getRealArray(itsAttr[MESHLENGTH]));
    ur->setMeshResolution(Attributes::getRealArray(itsAttr[MESHRESOLUTION]));
    ur->setTruncationOrder(Attributes::getReal(itsAttr[TRUNORDER]));
    ur->setSpaceCharge(Attributes::getReal(itsAttr[SPACECHARGE]));
    ur->setEmitParticles(Attributes::getReal(itsAttr[EMITPARTICLES]));
    ur->setTimeStepRatio(Attributes::getReal(itsAttr[TIMESTEPRATIO]));
    ur->setTotalTime(Attributes::getReal(itsAttr[TOTALTIME]));
    ur->setLFringe(Attributes::getReal(itsAttr[LFRINGE]));

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(ur);
}
