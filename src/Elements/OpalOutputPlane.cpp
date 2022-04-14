#include "Elements/OpalOutputPlane.h"
#include "Classic/AbsBeamline/OutputPlane.h"
#include "AbstractObjects/Attribute.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
//#include "Solvers/ScatteringPhysics.h"

std::string docstring = 
    std::string("The \"OUTPUTPLANE\" element writes out position at which")+
    std::string("trajectories cross a given plane.");

OpalOutputPlane::OpalOutputPlane():
    OpalElement(SIZE, "OUTPUTPLANE", docstring.c_str()) {

    itsAttr[CENTRE_X] = Attributes::makeReal
                      ("CENTRE_X", "x position of the plane centre [m]");
    itsAttr[CENTRE_Y] = Attributes::makeReal
                      ("CENTRE_Y", "y position of the plane centre [m]");
    itsAttr[CENTRE_Z] = Attributes::makeReal
                      ("CENTRE_Z", "z position of the plane centre [m]");
    itsAttr[NORMAL_X] = Attributes::makeReal
                      ("NORMAL_X", "x component of the normal to the plane");
    itsAttr[NORMAL_Y] = Attributes::makeReal
                      ("NORMAL_Y", "y component of the normal to the plane");
    itsAttr[NORMAL_Z] = Attributes::makeReal
                      ("NORMAL_Z", "z component of the normal to the plane");
    itsAttr[REFERENCE_ALIGNMENT_PARTICLE] = Attributes::makeReal
                      ("REFERENCE_ALIGNMENT_PARTICLE",
                       "Set to 1 to align output plane to reference particle");
    itsAttr[TOLERANCE] = Attributes::makeReal
                      ("TOLERANCE", "Tolerance on position of track intercept");
    itsAttr[WIDTH] = Attributes::makeReal
                      ("WIDTH", "Full width of the output plane, defined in the lab coordinate system (*not* the output plane coordinate system).");
    itsAttr[HEIGHT] = Attributes::makeReal
                      ("HEIGHT", "Full height of the output plane, defined in the lab coordinate system (*not* the output plane coordinate system)");
    itsAttr[RADIUS] = Attributes::makeReal
                      ("RADIUS", "Maximum distance from centre of plane for crossings");
    itsAttr[ALGORITHM] = Attributes::makeString
                      ("ALGORITHM",
                       "'EXTRAPOLATION' or 'RK4'");
    itsAttr[MATERIAL] = Attributes::makeString
                      ("MATERIAL", "Name of material that will be applied during crossing (default is no material)");
    itsAttr[MATERIAL_THICKNESS] = Attributes::makeReal
                      ("MATERIAL_THICKNESS", "Thickness of material that will be used during crossing (default is no material)");
    /*
    registerRealAttribute("CENTRE_X");
    registerRealAttribute("CENTRE_Y");
    registerRealAttribute("CENTRE_Z");
    registerRealAttribute("NORMAL_X");
    registerRealAttribute("NORMAL_Y");
    registerRealAttribute("NORMAL_Z");
    registerRealAttribute("REFERENCE_ALIGNMENT_PARTICLE");
    registerRealAttribute("TOLERANCE");
    registerRealAttribute("WIDTH");
    registerRealAttribute("HEIGHT");
    registerRealAttribute("RADIUS");
    registerStringAttribute("ALGORITHM");
    registerStringAttribute("MATERIAL");
    registerRealAttribute("MATERIAL_THICKNESS");
    */
    registerOwnership();

    //setElement((new OutputPlane("OUTPUTPLANE"))->makeAlignWrapper());
    setElement(new OutputPlane("OUTPUTPLANE"));
}


OpalOutputPlane::OpalOutputPlane(const std::string &name, OpalOutputPlane *parent):
    OpalElement(name, parent) {
    //setElement((new OutputPlane(name))->makeAlignWrapper());
    setElement(new OutputPlane(name));
}


OpalOutputPlane::~OpalOutputPlane() {
}


OpalOutputPlane *OpalOutputPlane::clone(const std::string &name) {
    return new OpalOutputPlane(name, this);
}

/*
void OpalOutputPlane::fillRegisteredAttributes(const ElementBase &base, ValueFlag flag) {
    OpalElement::fillRegisteredAttributes(base, flag);
}
*/


void OpalOutputPlane::update() {
    OpalElement::update();
    Vector_t centre(
        Attributes::getReal(itsAttr[CENTRE_X]),
        Attributes::getReal(itsAttr[CENTRE_Y]),
        Attributes::getReal(itsAttr[CENTRE_Z])
    );
    Vector_t normal(
        Attributes::getReal(itsAttr[NORMAL_X]),
        Attributes::getReal(itsAttr[NORMAL_Y]),
        Attributes::getReal(itsAttr[NORMAL_Z])
    );
    OutputPlane *output =
        //dynamic_cast<OutputPlane *>(getElement()->removeWrappers());
        dynamic_cast<OutputPlane *>(getElement());
    output->setCentre(centre);
    output->setNormal(normal);

    double tolerance = 1e-6;
    if (itsAttr[TOLERANCE]) {
        tolerance = Attributes::getReal(itsAttr[TOLERANCE]);
    }
    output->setTolerance(tolerance);

    if (itsAttr[WIDTH]) {
        double width = Attributes::getReal(itsAttr[WIDTH]);
        output->setHorizontalExtent(width/2.0);
    }
    if (itsAttr[HEIGHT]) {
        double height = Attributes::getReal(itsAttr[HEIGHT]);
        output->setVerticalExtent(height/2.0);
    }
    if (itsAttr[RADIUS]) {
        double radius = Attributes::getReal(itsAttr[RADIUS]);
        output->setRadialExtent(radius);
    }

    if (itsAttr[REFERENCE_ALIGNMENT_PARTICLE]) {
        double ref = Attributes::getReal(itsAttr[REFERENCE_ALIGNMENT_PARTICLE]);
        int refAlign = floor(ref + 0.5);
        output->setRecentre(refAlign);
    }
    std::string algorithm = "RK4";
    if (itsAttr[ALGORITHM]) {
        algorithm = Attributes::getString(itsAttr[ALGORITHM]);
    }
    if (algorithm == "RK4") {
        output->setAlgorithm(OutputPlane::RK4);
    } else if (algorithm == "INTERPOLATION") {
        output->setAlgorithm(OutputPlane::INTERPOLATION);
    } else {
        throw OpalException("OpalOutputPlane::update()",
                            "Did not recognise ALGORITHM '"+algorithm+"'");
    }
    /*
    if (itsAttr[MATERIAL]) {
        std::string material = Attributes::getString(itsAttr[MATERIAL]);
        double thickness = Attributes::getReal(itsAttr[MATERIAL_THICKNESS]);
        std::shared_ptr<CollimatorPhysics> physics;
        physics.reset(new CollimatorPhysics("Custom", output, material));
        output->setMaterial(physics, thickness);
    }
    */

    // Transmit "unknown" attributes.
    OpalElement::updateUnknown(output);
}