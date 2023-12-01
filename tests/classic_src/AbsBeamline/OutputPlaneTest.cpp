// Copyright (c) 2023, Chris Rogers
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include <cmath>
#include <fstream>
#include <sstream>

#include "gtest/gtest.h"
#include "opal_test_utilities/SilenceTest.h"

#include "Classic/Physics/Physics.h"
#include "Classic/AbsBeamline/OutputPlane.h"
#include "Classic/AbsBeamline/Component.h"


class MockDipole : public Component {
public:
    MockDipole(const std::string /*&name*/) : bfield_m(0, 1, 0), efield_m() {}
    NullField &getField() {return nullfield_m;}
    const NullField &getField() const {return nullfield_m;}
    void accept(BeamlineVisitor& /*&visitor*/) const {}
    void initialise(PartBunchBase<double, 3>* /*bunch*/, double& /*start*/, double& /*end*/) {}
    void finalise() {}
    bool bends() const {return true;}
    void getDimensions(double& /*zBegin*/, double& /*zEnd*/) const {}
    StraightGeometry& getGeometry() {return geom_m;}
    const StraightGeometry& getGeometry() const {return geom_m;}
    ElementBase* clone() const {return NULL;}

    void setField(Vector_t bfield, Vector_t efield) {bfield_m = bfield; efield_m = efield;}

    bool apply(const Vector_t& /*R*/, const Vector_t& /*P*/, const double& /*t*/, Vector_t& E, Vector_t& B) override {
        B = bfield_m;
        E = efield_m;
        return false;
    }
private:
    NullField nullfield_m;
    StraightGeometry geom_m;
    Vector_t bfield_m;
    Vector_t efield_m;
};

// Set also - lots of testing done in pyopal level
class OutputPlaneTest : public ::testing::Test {
public:
    OutputPlaneTest() : dipole_m("") {
    }

    void SetUp( ) {
        output_m.reset(new OutputPlane());
        output_m->setRecentre(-1);

    }

    void TearDown( ) {
    }

    ~OutputPlaneTest() {
        output_m.release();
    }

    std::unique_ptr<OutputPlane> output_m;
    MockDipole dipole_m;
    const double q2m = 1/Physics::c/Physics::c/Physics::m_p;

    // OpalTestUtilities::SilenceTest silencer_m;
};

TEST_F(OutputPlaneTest, TestSetGet) {
    // test basic set and get methods
    output_m->setGlobalFieldMap(&dipole_m);
    EXPECT_EQ(output_m->getGlobalFieldMap(), &dipole_m);

    output_m->setNormal(Vector_t(1.0, 2.0, 3.0));
    Vector_t ref = Vector_t(1.0, 2.0, 3.0)/std::sqrt(14);
    EXPECT_NEAR(output_m->getNormal()[0], ref[0], 1e-9);
    EXPECT_NEAR(output_m->getNormal()[1], ref[1], 1e-9);
    EXPECT_NEAR(output_m->getNormal()[2], ref[2], 1e-9);

    output_m->setCentre(Vector_t(4.0, 5.0, 6.0));
    EXPECT_EQ(output_m->getCentre(), Vector_t(4.0, 5.0, 6.0));

    output_m->setTolerance(1e-6);
    EXPECT_EQ(output_m->getTolerance(), 1e-6);
}


TEST_F(OutputPlaneTest, TestCheckOne_FieldOff) {
    //checkOne(const double tstep, const double chargeToMass,
    //         double& t, Vector_t& R, Vector_t& P)
    output_m->setGlobalFieldMap(&dipole_m);
    output_m->setCentre(Vector_t(0.0, 0.0, 1.0));
    output_m->setNormal(Vector_t(1.0, 0.0, 1.0));
    output_m->setTolerance(1e-9);

    Vector_t R, P;
    bool crossing;
    double t = 0.0;

    R = Vector_t(1.0, 0.0, 1.0);
    P = Vector_t(0.1, 0.0, 0.0);
    // c = 0.3 m/ns and particle is 1.0 m away; so should fail without doing RK4
    // step at all
    crossing = output_m->checkOne(0, 1.0e-9, 0.0, t, R, P);
    EXPECT_EQ(R, Vector_t(1.0, 0.0, 1.0));
    EXPECT_EQ(P, Vector_t(0.1, 0.0, 0.0));
    EXPECT_FALSE(crossing);
    // betagamma = 0.1 and particle is ~0.07 m away; so should fail after trying
    // a single RK4 step (which goes about 0.03 m)
    R = Vector_t(1-0.05, 0.0, -0.05);
    crossing = output_m->checkOne(0, 1.0e-9, 0.0, t, R, P);
    EXPECT_EQ(R, Vector_t(1.0-0.05, 0.0, -0.05));
    EXPECT_EQ(P, Vector_t(0.1, 0.0, 0.0));
    EXPECT_FALSE(crossing);
    // betagamma = 0.1 and particle is 0.014 m away; so should succeed after a
    // couple of RK4 step (which goes about 0.03 m)
    R = Vector_t(1.0-0.01, 0.0, -0.01);
    crossing = output_m->checkOne(0, 1.0, 0.0, t, R, P);
    EXPECT_NEAR(R[0], 1.01, 1e-5);
    EXPECT_NEAR(R[1], 0.0, 1e-5);
    EXPECT_NEAR(R[2], -0.01, 1e-5);
    EXPECT_EQ(P, Vector_t(0.1, 0.0, 0.0));
    EXPECT_TRUE(crossing);
    // particle is past the plane but close; should fail after trying a single
    // RK4 step
    R = Vector_t(1.0+0.01, 0.0, +0.01);
    crossing = output_m->checkOne(0, 1.0, 0.0, t, R, P);
    EXPECT_NEAR(R[0], 1.0+0.01, 1e-5);
    EXPECT_NEAR(R[1], 0.0, 1e-5);
    EXPECT_NEAR(R[2], +0.01, 1e-5);
    EXPECT_EQ(P, Vector_t(0.1, 0.0, 0.0));
    EXPECT_FALSE(crossing);
    // particle is past the plane but travelling backwards; should pass
    R = Vector_t(1.0+0.01, 0.0, +0.01);
    P = Vector_t(-0.1, 0.0, 0.0);
    crossing = output_m->checkOne(0, 1.0, 0.0, t, R, P);
    EXPECT_NEAR(R[0], 1.0-0.01, 1e-5);
    EXPECT_NEAR(R[1], 0.0, 1e-5);
    EXPECT_NEAR(R[2], +0.01, 1e-5);
    EXPECT_EQ(P, Vector_t(-0.1, 0.0, 0.0));
    EXPECT_TRUE(crossing);
}

double magnitude(Vector_t v) {
    return std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}


// ps: trying to figure out OPAL units, I put some print statements in ParallelCyclotronTracker
//
// itsBunch_m->Q[i]/itsBunch_m->M[i]*Physics::c*Physics::c/Physics::q_e/1e3 = 9.57883e+13
// std::cerr << "CHARGE TO MASS " << i << " Qi " << itsBunch_m->Q[i] << " Mi "
//           << itsBunch_m->M[i] << " c " << Physics::c << " qe " << Physics::q_e << std::endl;
// CHARGE TO MASS 0 Qi 1.60218e-19 Mi 0.938272 c 2.99792e+08 qe 1.60218e-19

TEST_F(OutputPlaneTest, TestRK4StepNoField) {
	Vector_t pVecRef(-1.0, -0.0, -0.4); // moving at beta gamma = 1.414
	Vector_t rVecRef(2., 3., 4.);
    Vector_t rVec = rVecRef;
    Vector_t pVec = pVecRef;
	double time = 2.0e-9;
    double tstep = 1.0e-9;
    double gammabeta = magnitude(pVecRef);
    double beta = gammabeta/std::sqrt(1+gammabeta*gammabeta);
    double gamma = gammabeta/beta;


    // field not set
    EXPECT_THROW(output_m->RK4Step(tstep, q2m, time, rVec, pVec), GeneralClassicException);

    // zero field
    Vector_t speed = pVecRef/gamma*Physics::c; // (beta gamma)* c / gamma = velocity
    rVecRef = speed*tstep+rVec;
    std::shared_ptr<MockDipole> testDipole;
    testDipole.reset(new MockDipole("mock"));
    testDipole->setField(Vector_t(), Vector_t());

    output_m->setGlobalFieldMap(testDipole.get());
    output_m->RK4Step(tstep, q2m, time, rVec, pVec);
    EXPECT_EQ(pVec, pVecRef);
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(rVec[i], rVecRef[i], 1e-6);
    }
}

TEST_F(OutputPlaneTest, TestRK4StepDipole) {
    // dipole in z direction
    // B = 2.0 Tesla
    // p = 1.0 GeV/c
    // rho = 1e9/2/c
    double time = 0.0;
    double tstep = 1.0e-2;

    std::shared_ptr<MockDipole> testDipole;
    testDipole.reset(new MockDipole("mock"));
    testDipole->setField(Vector_t(0.0, 0.0, -2.0), Vector_t());
    Vector_t rVecRef = Vector_t(0.0, 0.0, 0.0);
    Vector_t pVecRef = Vector_t(0.0, 1.0/Physics::m_p, 1.0);
    Vector_t rVec = rVecRef;
    Vector_t pVec = pVecRef;
    output_m->setGlobalFieldMap(testDipole.get());
    output_m->RK4Step(tstep, q2m, time, rVec, pVec);
    EXPECT_NEAR(pVecRef[2], pVec[2], 1e-9);
    EXPECT_NEAR(pVecRef[0]*pVecRef[0]+pVecRef[1]*pVecRef[1],
                pVec[0]*pVec[0]+pVec[1]*pVec[1], 1e-8);
}

TEST_F(OutputPlaneTest, TestRK4StepEField) {
    // dipole in z direction
    // B = 2.0 Tesla
    // p = 1.0 GeV/c
    // rho = 1e9/2/c
    double time = 0.0;
    double tstep = 1.0e-9;

    std::shared_ptr<MockDipole> testDipole;
    testDipole.reset(new MockDipole("mock"));
    testDipole->setField(Vector_t(), Vector_t(0.0, 0.0, 2.0));
    Vector_t rVecRef = Vector_t(1e9/2/Physics::c, 0.0, 0.0);
    Vector_t pVecRef = Vector_t(0.0, 1.0/0.938272, 0.0);
    Vector_t rVec = rVecRef;
    Vector_t pVec = pVecRef;
    output_m->setGlobalFieldMap(testDipole.get());
    output_m->RK4Step(tstep, q2m, time, rVec, pVec);
    EXPECT_NEAR(pVec[0], pVecRef[0], 1e-9);
    EXPECT_NEAR(pVec[1], pVecRef[1], 1e-9);
}

TEST_F(OutputPlaneTest, TestCheckOne_FieldOn) {
    //checkOne(const double tstep, const double chargeToMass,
    //         double& t, Vector_t& R, Vector_t& P)
    output_m->setGlobalFieldMap(&dipole_m);
    output_m->setCentre(Vector_t(0.0, 0.0, 0.0));
    output_m->setNormal(Vector_t(1.0, 0.0, 0.0));
    output_m->setTolerance(1e-9);
    output_m->setAlgorithm(OutputPlane::algorithm::RK4STEP);

    Vector_t R, P;
    bool crossing;
    double q2m=1.0*Physics::q_e; // proton mass is 0.938
    // mass excluded since RK4 solvers works in terms of momentum
    double t = 0.0;

    // betagamma = 0.1 and particle is 0.02 m away; so should succeed after a
    // couple of RK4 step (which goes about 0.03 m)
    R = Vector_t(-0.02, 0.0, 0.0);
    P = Vector_t(0.1, 0.1, 0.0);
    crossing = output_m->checkOne(0, 1.0, q2m, t, R, P);
    
    EXPECT_NEAR(R[0], 0.0, 1e-5);
    EXPECT_NEAR(R[1], 0.02, 1e-5);
    EXPECT_NEAR(R[2], 0.0, 1e-5);
    EXPECT_NEAR(P[0], 0.1, 1e-5);
    EXPECT_NEAR(P[1], 0.1, 1e-5);
    EXPECT_NEAR(P[2], 0.0, 1e-5);
    
    EXPECT_TRUE(crossing);
}

TEST_F(OutputPlaneTest, TestCentre) {
    // set centre
    output_m->setCentre(Vector_t(2.0, 3.0, 4.0));
    output_m->setNormal(Vector_t(1.0, 0.0, 0.0));
    output_m->setTolerance(1e-9);

    Vector_t R, P, normP, R1, R2, a, b, normaltoplane, unitnormaltoplane;
    R = Vector_t(-0.02, 5.0, 7.1);
    P = Vector_t(0.5, 0.8, 0.3);

    //re centre plane
    output_m->recentre(R, P);

    double Ptot = std::sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
    normP = Vector_t(P[0]/Ptot, P[1]/Ptot, P[2]/Ptot);

    // check plane is recentred properly
    EXPECT_EQ(output_m->getCentre(), R);
    EXPECT_EQ(output_m->getNormal(), normP);
 }

TEST_F(OutputPlaneTest, TestExtent) {
    output_m->setCentre(Vector_t(2.0, 3.0, 4.0));
    output_m->setNormal(Vector_t(1.0, 0.0, 0.0));
    output_m->setTolerance(1e-9);
    double halfwidth, halfheight;

    halfheight = 10;
    halfwidth = 10;

    double horizextent = output_m->getHorizontalExtent();
    // calling getHorizontalExtent before settingthe horizontal extent returns a value of -1
    EXPECT_EQ(-1, horizextent);

    output_m->setHorizontalExtent(halfwidth);
    output_m->setVerticalExtent(halfheight);
    double vertextent = output_m->getVerticalExtent();
    horizextent = output_m->getHorizontalExtent();

    EXPECT_EQ( horizextent, halfwidth);
    EXPECT_EQ( vertextent, halfheight);

}


TEST_F(OutputPlaneTest, checkOneTest){
    OutputPlane* myOutputPlane = new OutputPlane("TestPlane"); 
    MockDipole* myDipole = new MockDipole("TestDipole");
    // check that checkOne returns true when the particle is within one time step
    // of the OutputPlane and false when it is not within one timestep
    myOutputPlane->setGlobalFieldMap(myDipole);

    Vector_t R, P;
    Vector_t R1, P1;

    double tStep = 1.0e-9; // s

    Vector_t normal = {1.0, 0.0, 1.0};
    Vector_t centre = {0.0, 0.0, 1.0};
    myOutputPlane->setCentre(centre);
    myOutputPlane->setNormal(normal);
    myOutputPlane->setTolerance(1e-9);

    R = Vector_t(0.1, 0.0, 0.1);
    P = Vector_t(0.1, 0.0, 0.0); 
    double t = 0.0;
    bool cross;
    cross = myOutputPlane->checkOne(0, tStep, 0.0, t, R, P);
    EXPECT_FALSE(cross);

    // opposite case where it is near to the plane
    // distance from plane
    R1 = Vector_t(1.0-0.01, 0.0, -0.01);
    P1 = Vector_t(0.1, 0.0, 0.0); 

    bool cross1;
    cross1 = myOutputPlane->checkOne(0, tStep, 0.0, t, R1, P1);
    EXPECT_TRUE(cross1);
    delete myDipole;
    delete myOutputPlane;
}
