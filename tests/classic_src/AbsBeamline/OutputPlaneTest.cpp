#include <cmath>
#include <fstream>
#include <sstream>

#include "gtest/gtest.h"
#include "opal_test_utilities/SilenceTest.h"
#include "classic_src/AbsBeamline/TrackingTest/Tracking.h"

#include "Classic/Physics/Physics.h"
#include "Classic/AbsBeamline/OutputPlane.h"
#include "Classic/AbsBeamline/Component.h"


class MockDipole : public Component {
public:
    MockDipole(const std::string /*&name*/) {}
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

    bool apply(const Vector_t& /*R*/, const Vector_t& /*P*/, const double& /*t*/, Vector_t& E, Vector_t& B) override {
        B[0] = 0.;
        B[1] = 1.;
        B[2] = 0.;        
        
        E[0] = 0.;
        E[1] = 0.;
        E[2] = 0.;
        return false;
    }
private:
    NullField nullfield_m;
    StraightGeometry geom_m;
};


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
	const double m_p=938.2720813; // proton mass MeV/c^2

    // OpalTestUtilities::SilenceTest silencer_m;
};

TEST_F(OutputPlaneTest, TestSetGet) {
    // test basic set and get methods
    output_m->setGlobalFieldMap(&dipole_m);
    EXPECT_EQ(output_m->getGlobalFieldMap(), &dipole_m);

    output_m->setNormal(Vector_t(1.0, 2.0, 3.0));
    Vector_t ref = Vector_t(1.0, 2.0, 3.0)/sqrt(14);
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
    crossing = output_m->checkOne(0, 1.0, 0.0, t, R, P);
    EXPECT_EQ(R, Vector_t(1.0, 0.0, 1.0));
    EXPECT_EQ(P, Vector_t(0.1, 0.0, 0.0));
    EXPECT_FALSE(crossing);
    // betagamma = 0.1 and particle is ~0.07 m away; so should fail after trying
    // a single RK4 step (which goes about 0.03 m)
    R = Vector_t(1-0.05, 0.0, -0.05);
    crossing = output_m->checkOne(0, 1.0, 0.0, t, R, P);
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

TEST_F(OutputPlaneTest, TestRK4Step) {
	Vector_t pvec(1000., 0., 0.);
	Vector_t rvec(0., 0., 0.);
	double time = 0.;

	double p = euclidean_norm(pvec);
	double energy = ::sqrt(p*p+m_p*m_p);

	double x_in[] = {time, rvec[0], rvec[1], rvec[2], energy, pvec[0], pvec[1], pvec[2]};
	Tracking tracking;
	tracking.setCharge(1.);
	tracking.setMaximumStepNumber(1);
	for (size_t i = 0; i < 8; ++i) {
		std::cerr << x_in[i] << "   ";
	}
	std::cerr << std::endl;
	try {
		tracking.integrate(Tracking::t, 1., x_in, &dipole_m, 1); // 1 ns step
	} catch (GeneralClassicException exc) {
		std::cerr << exc.what() << std::endl;
	}
	for (size_t i = 0; i < 8; ++i) {
		std::cerr << x_in[i] << "   ";
	}
	std::cerr << std::endl;

	pvec /= m_p;

	std::cerr << time << "   ";
	for (size_t i = 0; i < 3; ++i) {
		std::cerr << rvec[i] << "   ";
	}
	std::cerr << ::sqrt(dot(pvec, pvec)+1)*m_p << "   ";
	for (size_t i = 0; i < 3; ++i) {
		std::cerr << m_p*pvec[i] << "   ";
	}
	std::cerr << std::endl;

	output_m->setGlobalFieldMap(&dipole_m);
    output_m->RK4Step(1., 1./m_p, time, rvec, pvec);
	std::cerr << time << "   ";
	for (size_t i = 0; i < 3; ++i) {
		std::cerr << rvec[i] << "   ";
	}
	std::cerr << ::sqrt(dot(pvec, pvec)+1)*m_p << "   ";
	for (size_t i = 0; i < 3; ++i) {
		std::cerr << m_p*pvec[i] << "   ";
	}
	std::cerr << std::endl;
}


TEST_F(OutputPlaneTest, TestCheckOne_FieldOn) {
    //checkOne(const double tstep, const double chargeToMass,
    //         double& t, Vector_t& R, Vector_t& P)
    output_m->setGlobalFieldMap(&dipole_m);
    output_m->setCentre(Vector_t(0.0, 0.0, 0.0));
    output_m->setNormal(Vector_t(1.0, 0.0, 0.0));
    output_m->setTolerance(1e-9);
    output_m->setAlgorithm(OutputPlane::algorithm::RK4);

    Vector_t R, P;
    bool crossing;
    double q2m=1.0*Physics::q_e; // proton mass is 0.938
    // mass excluded since RK4 solvers works in terms of momentum
    double t = 0.0;

    /*
    R = Vector_t(-1.0, 0.0, 0.0);
    P = Vector_t(0.1, 0.0, 0.0);
    // betagamma = 0.1 and particle is 0.05 m away; so should fail after trying
    // a single RK4 step (which goes about 0.03 m)
    R = Vector_t(-0.05, 0.0, 0.0);
    crossing = output_m->checkOne(0, 1.0, q2m, t, R, P);
    EXPECT_EQ(R, Vector_t(-0.05, 0.0, 0.0));
    EXPECT_EQ(P, Vector_t(0.1, 0.0, 0.0));
    EXPECT_FALSE(crossing);
    */

    // betagamma = 0.1 and particle is 0.02 m away; so should succeed after a
    // couple of RK4 step (which goes about 0.03 m)
    R = Vector_t(-0.02, 0.0, 0.0);
    P = Vector_t(0.1, 0.1, 0.0);
    crossing = output_m->checkOne(0, 1.0, q2m, t, R, P);
    
    EXPECT_NEAR(R[0], 0.0, 1e-5);
    EXPECT_NEAR(R[1], 0.02001200, 1e-5);
    EXPECT_NEAR(R[2], 0.00060012, 1e-5);
    EXPECT_NEAR(P[0], 0.0998201, 1e-5);
    EXPECT_NEAR(P[1], 0.1, 1e-5);
    EXPECT_NEAR(P[2], 0.00599585, 1e-5);    
    
    EXPECT_TRUE(crossing);
}

TEST_F(OutputPlaneTest, TestEField) {
    //checkOne(const double tstep, const double chargeToMass,
    //         double& t, Vector_t& R, Vector_t& P)


    output_m->setGlobalFieldMap(&dipole_m);
    output_m->setCentre(Vector_t(0.0, 0.0, 0.0));
    output_m->setNormal(Vector_t(1.0, 0.0, 0.0));
    output_m->setTolerance(1e-9);
    output_m->setAlgorithm(OutputPlane::algorithm::RK4);

    Vector_t R, P;
    double q2m=1.0*Physics::q_e; // proton mass is 0.938
    // mass excluded since RK4 solvers works in terms of momentum
    double t = 0.0;


    R = Vector_t(0.0, 0.0, 0.0); 
    P = Vector_t(0, 0.1, 0.0);
    output_m->RK4Step(1.0, q2m, t, R, P);

    // when E[0] in mock dipole is set to 10 the follwing checks should pass
    //EXPECT_NEAR(P[0], 2.99792, 1e-5);
    //EXPECT_NEAR(P[1], 0.1, 1e-5);
    //EXPECT_NEAR(P[2], 0.0, 1e-5); 

    //EXPECT_NEAR(R[0], 0.21709, 1e-5);
    //EXPECT_NEAR(R[1], 0.03177, 1e-5);
    //EXPECT_NEAR(R[2], 0.0, 1e-5); 
    
    // when E[0] is set to zero there should be no change in P[0] 
    EXPECT_NEAR(P[0], 0.0, 1e-5);
    EXPECT_NEAR(P[1], 0.1, 1e-5);
    EXPECT_NEAR(P[2], 0.0, 1e-5); 

    
    EXPECT_NEAR(R[0], 0.0, 1e-5);
    EXPECT_NEAR(R[1], 0.03177, 1e-5);
    EXPECT_NEAR(R[2], 0.0, 1e-5); 

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

    double Ptot = sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
    normP = Vector_t(P[0]/Ptot, P[1]/Ptot, P[2]/Ptot);

    // check plane is recentred properly
    EXPECT_EQ(output_m->getCentre(), R);
    EXPECT_EQ(output_m->getNormal(), normP);
    
    /*
    // check that plane defined by R and another 2 points on the plane is parallel to P(ie normal to the plane)
    R1 = Vector_t((R[0]-1), (R[1]-2), (R[2]-3));
    R2 = Vector_t((R[0]-4), (R[1]-5), (R[2]-6));

    a = R1 - R;
    b = R2 - R;

    normaltoplane = Vector_t(a[1]*b[2] - b[1]*a[2], -(a[0]*b[2] - b[0]*a[2]), a[0]*b[1] - b[0]*a[1]);

    double normaltoplanetot = sqrt(normaltoplane[0]*normaltoplane[0] + normaltoplane[1]*normaltoplane[1] + normaltoplane[2]*normaltoplane[2]);
    unitnormaltoplane = Vector_t(normaltoplane[0]/normaltoplanetot, normaltoplane[1]/normaltoplanetot, normaltoplane[2]/normaltoplanetot);

    EXPECT_NEAR(1.0, unitnormaltoplane[0]*normP[0] + unitnormaltoplane[1]*normP[1] + unitnormaltoplane[2]*normP[2], 1e-5);
    */
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

    double tStep = 1.0; // ns

    Vector_t normal = {1.0, 0.0, 1.0};
    Vector_t centre = {0.0, 0.0, 1.0};
    myOutputPlane->setCentre(centre);
    myOutputPlane->setNormal(normal);
    myOutputPlane->setTolerance(1e-9);

    R = Vector_t(0.1, 0.0, 0.1);
    P = Vector_t(0.1, 0.0, 0.0); 
    double t = 0.0;
    bool cross;
    // distance from plane

    //Vector_t delta = R - centre;
    
    cross = myOutputPlane->checkOne(0, tStep, 0.0, t, R, P);

    //double dist = normal[0]*delta[0] + normal[1]*delta[1] + normal[2]*delta[2];
    //double sStep = tStep*euclidean_norm(P)*Physics::c*1e-9;

    //myOutputPlane->checkOne(0, 1.0, 1.1, t, R, P);
    //bool val;
    //if (fabs(dist) > sStep) {
    //    val = false;
    //} else {
    //    val = true;
    //}
    EXPECT_FALSE(cross);

    // opposite case where it is near to the plane
    //R1 = Vector_t(0.01, 0.01, 0.02);
    //P1 = Vector_t(0,0, 0); //  initial momnetum
    //double t = 10.0;

    // distance from plane
    R1 = Vector_t(1.0-0.01, 0.0, -0.01);
    P1 = Vector_t(0.1, 0.0, 0.0); 

    bool cross1;

    //Vector_t delta1 = R1 - centre;

    //double dist1 = normal[0]*delta[0] + normal[1]*delta[1] + normal[2]*delta[2];
    //double sStep1 = tStep1*euclidean_norm(P)*Physics::c*1e-9;
    cross1 = myOutputPlane->checkOne(0, tStep, 0.0, t, R1, P1);
    //EXPECT_TRUE(false) << cross1;
    //myOutputPlane->checkOne(0, 1.0, 1.1, t, R, P);

    //if (fabs(dist1) > sStep1) {
    //    val = false;
    //}else{
    //    val = true;
    //}
    EXPECT_TRUE(cross1);

}
