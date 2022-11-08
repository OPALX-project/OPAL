#include <sstream>

#include "gtest/gtest.h"
#include "Attributes/Attributes.h"
#include "Classic/AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "Classic/AbsBeamline/EndFieldModel/Tanh.h"
#include "Classic/AbsBeamline/ScalingFFAMagnet.h"
#include "Elements/OpalScalingFFAMagnet.h"
#include "Physics/Units.h"
#include "Physics/Physics.h"

#include "opal_test_utilities/SilenceTest.h"

void setReal(OpalScalingFFAMagnet& mag, std::string attName, double value) {
    Attribute* att = mag.findAttribute(attName);
    ASSERT_NE(att, nullptr);
    Attributes::setReal(*att, value);
}

TEST(OpalScalingFFAMagnetTest, TestConstructorDestructor) {
    OpalTestUtilities::SilenceTest silencer;

    OpalScalingFFAMagnet mag1;
    EXPECT_EQ(mag1.getOpalName(), "SCALINGFFAMAGNET");
}

TEST(OpalScalingFFAMagnetTest, TestUpdate) {
    OpalScalingFFAMagnet opalMag1;
    std::cerr << "HELLO 0" << std::endl;
    setReal(opalMag1, "B0", 1.0);
    setReal(opalMag1, "R0", 12.5);
    setReal(opalMag1, "FIELD_INDEX", 2.0);
    setReal(opalMag1, "TAN_DELTA", 0.7);
    setReal(opalMag1, "MAX_Y_POWER", 3.0);
    setReal(opalMag1, "END_LENGTH", 0.8);
    setReal(opalMag1, "CENTRE_LENGTH", 5.0);
    setReal(opalMag1, "RADIAL_NEG_EXTENT", 0.55);
    setReal(opalMag1, "RADIAL_POS_EXTENT", 0.65);
    setReal(opalMag1, "HEIGHT", 6.0);
    std::cerr << "HELLO 0a" << std::endl;
    opalMag1.update();
    std::cerr << "HELLO 0b" << std::endl;
    ElementBase* element = opalMag1.getElement();
    ScalingFFAMagnet* mag1 = dynamic_cast<ScalingFFAMagnet*>(element);
    std::cerr << "HELLO 0c" << std::endl;
    mag1->setupEndField();
    ASSERT_NE(mag1, nullptr);
    std::cerr << "HELLO 1" << std::endl;
    EXPECT_NEAR(mag1->getDipoleConstant(), 1.0*Units::T2kG, 1e-9);
    EXPECT_NEAR(mag1->getR0(), 12.5*Units::m2mm, 1e-9);
    EXPECT_NEAR(mag1->getFieldIndex(), 2.0, 1e-9);
    EXPECT_NEAR(mag1->getTanDelta(), 0.7, 1e-9);
    EXPECT_EQ(mag1->getMaxOrder(), 3u);
    EXPECT_NEAR(mag1->getRMin(), 11.95*Units::m2mm, 1e-9);
    EXPECT_NEAR(mag1->getRMax(), 13.15*Units::m2mm, 1e-9);
    EXPECT_NEAR(mag1->getVerticalExtent(), 3.0*Units::m2mm, 1e-9);
    std::cerr << "HELLO 2" << std::endl;

    endfieldmodel::EndFieldModel* model = mag1->getEndField();
    endfieldmodel::Tanh* tanh = dynamic_cast<endfieldmodel::Tanh*>(model);
    ASSERT_NE(tanh, nullptr);
    EXPECT_NEAR(tanh->getX0(), 5.0/12.5/2, 1e-9);
    EXPECT_NEAR(tanh->getLambda(), 0.8/12.5, 1e-9);
    EXPECT_NEAR(mag1->getPhiStart(), 0.5*(4*0.8+5.0)/12.5, 1e-9);
    EXPECT_NEAR(mag1->getPhiEnd(), (4*0.8+5.0)/12.5, 1e-9);
    EXPECT_NEAR(mag1->getAzimuthalExtent(), (5*0.8+5.0/2.0)/12.5, 1e-9);
    std::cerr << "HELLO 3" << std::endl;

    setReal(opalMag1, "MAGNET_START", 7.0);
    setReal(opalMag1, "MAGNET_END", 8.0);
    setReal(opalMag1, "AZIMUTHAL_EXTENT", 9.0);
    opalMag1.update();
    mag1->setupEndField();
    EXPECT_NEAR(mag1->getPhiStart(), (7.0+5.0/2)/12.5, 1e-9);
    EXPECT_NEAR(mag1->getPhiEnd(), 8.0/12.5, 1e-9);
    EXPECT_NEAR(mag1->getAzimuthalExtent(), 9.0/12.5, 1e-9);
}

// r in mm, phi in rad; note the weird coordinate system
Vector_t cartesianCoord(double r, double phi) {
    phi -= Physics::pi/2.0;
    return Vector_t(-r*sin(phi), 0.0, r*cos(phi))-Vector_t(r, 0, 0);
}

TEST(OpalScalingFFAMagnetTest, TestFieldCheck) {
    OpalScalingFFAMagnet opalMag1;
    double r0 = 4.0;
    double b0 = 1.1;
    double PI = Physics::pi;
    setReal(opalMag1, "B0", b0);
    setReal(opalMag1, "R0", r0);
    setReal(opalMag1, "FIELD_INDEX", 2.0);
    setReal(opalMag1, "MAX_Y_POWER", 3.0);
    setReal(opalMag1, "END_LENGTH", r0*PI/64.0);
    setReal(opalMag1, "CENTRE_LENGTH", 16*r0*PI/64.0); // 1/8 of a circle
    setReal(opalMag1, "RADIAL_NEG_EXTENT", 2.0);
    setReal(opalMag1, "RADIAL_POS_EXTENT", 2.0);
    setReal(opalMag1, "HEIGHT", 6.0);
    setReal(opalMag1, "MAGNET_START", r0*12*PI/64.0);
    setReal(opalMag1, "MAGNET_END", 32*r0*PI/64.0);
    setReal(opalMag1, "AZIMUTHAL_EXTENT", r0*(16.0*PI/64.0));
    opalMag1.update();
    ElementBase* element = opalMag1.getElement();
    ScalingFFAMagnet* mag1 = dynamic_cast<ScalingFFAMagnet*>(element);
    mag1->setupEndField();
    ASSERT_NE(mag1, nullptr);
    // I just want to check that the bounding boxes/etc make sense
    double b0kG = b0*Units::T2kG;
    std::vector<double> position = {3.99, 4.01, 12, 20, 28, 35.99, 36.01};
    std::vector<bool> oob = {1, 0, 0, 0, 0, 0, 1};
    std::vector<double> field = {0.0, 0.0, b0kG/2, b0kG, b0kG/2, 0.0, 0.0};
    for (size_t i = 0; i < position.size(); ++i) {
        double r0mm = r0*Units::m2mm;
        double phi = position[i]*PI/64.0;
        Vector_t B, P, E;
        double t = 0.;
        // magnet start plus half centre length should have maximum field
        Vector_t rMiddle = cartesianCoord(r0mm, phi);
        bool outOfBounds = mag1->apply(rMiddle, P, t, E, B);
        EXPECT_EQ(outOfBounds, oob[i]);
        EXPECT_NEAR(B[1], field[i], 1e-4) << "failed for phi "
                                        << position[i] << "*PI/64 " << std::endl;
    }

    Euclid3D delta = mag1->getGeometry().getTotalTransform();
    Vector3D vec = delta.getVector();
    Vector3D rot = delta.getRotation().getAxis();
    EXPECT_EQ(vec(0), r0*Units::m2mm*(cos(mag1->getPhiEnd())-1));
    EXPECT_EQ(vec(1), 0.);
    EXPECT_EQ(vec(2), r0*Units::m2mm*sin(mag1->getPhiEnd()));

    EXPECT_EQ(rot(0), 0.);
    EXPECT_EQ(rot(1), -mag1->getPhiEnd());
    EXPECT_EQ(rot(2), 0.);

}