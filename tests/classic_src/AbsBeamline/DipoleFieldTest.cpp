#include "gtest/gtest.h"
#include <gtest/gtest_prod.h>
#include "AbsBeamline/Bend2D.h"
#include "Fields/Fieldmap.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "BeamlineCore/RBendRep.h"
#include "BeamlineCore/SBendRep.h"
#include "Algorithms/PartBunch.h"
#include "BeamlineCore/MultipoleRep.h"
#include "AbsBeamline/ElementBase.h"

#include "opal_test_utilities/SilenceTest.h"

#include <fstream>
using namespace std;

vector< vector<double> > partialsDerivB(const Vector_t &R,const Vector_t /*B*/, double stepSize, SBendRep* dummyField)
{
    // builds a matrix of all partial derivatives of B -> dx_i B_j
    vector< vector<double> > allPartials(3, vector<double>(3));
    double t = 0 ;
    Vector_t P, E;
    for(int i = 0; i < 3; i++)
	{
	  // B at the previous and next grid points R_prev,  R_next
	  Vector_t R_prev = R, R_next = R;
	  R_prev[i] -= stepSize;
	  R_next[i] += stepSize;
	  Vector_t B_prev, B_next;
	  dummyField->apply(R_prev, P, t, E, B_prev);
	  dummyField->apply(R_next, P, t, E, B_next);
	  for(int j = 0; j < 3; j++)
	    allPartials[i][j] = (B_next[j] - B_prev[j]) / (2 * stepSize);
	}
     return allPartials;
}

vector< vector<double> > partialsDerivB_5(const Vector_t &R,const Vector_t /*B*/, double stepSize, SBendRep* dummyField)
{
    // builds a matrix of all partial derivatives of B -> dx_i B_j
    vector< vector<double> > allPartials(3, vector<double>(3));
    double t = 0 ;
    Vector_t P, E;
    for(int i = 0; i < 3; i++)
	{
	  // B at the previous and next grid points R_prev,  R_next
	  Vector_t R_pprev = R, R_prev = R, R_next = R, R_nnext = R;
	  R_pprev(i) -= 2 * stepSize;
	  R_nnext(i) += 2 * stepSize;
	  R_prev(i) -= stepSize;
	  R_next(i) += stepSize;
	  Vector_t B_prev, B_next, B_pprev, B_nnext;
	  dummyField->apply(R_prev, P, t, E, B_prev);
	  dummyField->apply(R_next, P, t, E, B_next);
	  dummyField->apply(R_pprev, P, t, E, B_pprev);
	  dummyField->apply(R_nnext, P, t, E, B_nnext);
	  for(int j = 0; j < 3; j++)
	    allPartials[i][j] = (B_pprev[j] - 8 * B_prev[j] + 8 * B_next[j] - B_nnext[j]) / (12 * stepSize);
	}
     return allPartials;
}

double calcDivB(Vector_t &R, Vector_t B, double stepSize, SBendRep* dummyField )
{
    double div = 0;
    vector< vector<double> > partials (3, vector<double>(3));
    partials = partialsDerivB(R, B, stepSize, dummyField);
    for(int i = 0; i < 3; i++)
        div += partials[i][i];
    return div;
}

vector<double> calcCurlB(Vector_t &R, Vector_t B, double stepSize, SBendRep* dummyField)
{
    vector<double> curl(3);
    vector< vector<double> > partials(3, vector<double>(3));
    partials = partialsDerivB(R, B, stepSize, dummyField);
    curl[0] = (partials[1][2] - partials[2][1]);
    curl[1] = (partials[2][0] - partials[0][2]);
    curl[2] = (partials[0][1] - partials[1][0]);
    return curl;
}


TEST(Maxwell, Zeros)
{
    OpalTestUtilities::SilenceTest silencer;

    SBendRep* myMagnet = new SBendRep("myMagnet");
    myMagnet->BendBase::setFieldMapFN("1DPROFILE1-DEFAULT");
    myMagnet->BendBase::setLength(0.2);
    myMagnet->BendBase::setDesignEnergy(10.0e6);
    myMagnet->BendBase::setBendAngle(0.523599);//30 degrees
    myMagnet->BendBase::setFullGap(0.04);
    myMagnet->setElementPosition(0.5);
    //myMagnet->BendBase::setFieldAmplitude(-0.0350195, 0);
    myMagnet->BendBase::setEntranceAngle(0.);
    myMagnet->setK1(0); //set quadrupole component
    //the following are used to initialise myMagnet
    PartData* partData = new PartData();
    PartBunch* bunch = new PartBunch(partData);
    bunch->resetM(0.938);
    bunch->setdT(1.0e-12);//time step
    double startField = 2.0, endField = 10.0 ;
    myMagnet->Bend2D::initialise(bunch, startField, endField);
    delete partData;
    delete bunch;
    Vector_t R(0. ,0., 0.);
    Vector_t P(0.), E(0.);
    //double radius = 1;
    double stepSize = 1.0e-6, x, z;
    //int step;
    int counter = 0;
    //ofstream fout("some_data");
    for(z = 0.0; z <0.0015; z+= 0.0015)
        for(x = 0.; x<0.04; x += 0.04)
	  for(double phi = -Physics::pi / 7.1 ; phi < 2/3. * Physics::pi; phi += Physics::pi/2000.)
                {
		  // step = phi/(Physics::pi/20);
		  //std::cout<<"Step #"<<step<<endl;
		  counter ++;
		  Vector_t B(0.0);
		  R(0) = (myMagnet->Bend2D::designRadius_m + x) * cos(phi);
		  R(1) = z;
		  R(2) = (myMagnet->Bend2D::designRadius_m + x) * sin(phi);
		  double t = 0;
		  myMagnet->apply(R, P, t , E, B);
		  //B /= myMagnet->fieldAmplitude_m; //normalisation
		  //fout<<phi<<' '<<B[1] / myMagnet->fieldAmplitude_m<<endl;
		  //myMagnet.Bend::calculateMapField(R, B);
		  //std::cout<< "Position: " <<"phi="<<phi<<" x="<<x<<" z="<<z<<endl;
		  //std::cout<< "Field:" <<' '<<B[0]<<' ' <<B[1]<<' '<<B[2]<<endl;
		  double div = 0;
		  div = calcDivB(R, B, stepSize, myMagnet);
		  //fout<<phi<<' '<<z<<' '<<div<<' '<<endl;
		  vector<double> curl;
		  EXPECT_NEAR(div, 0.0, 0.15);
		  curl = calcCurlB(R, B, stepSize, myMagnet);
		  for(int k=0; k<3; k++) curl[k] /= myMagnet->fieldAmplitude_m;
		  //fout<<phi<<' '<<z<<' '<<sqrt(pow(curl[0], 2) + pow(curl[1], 2) + pow(curl[2], 2))<<endl;
		  //fout<<phi<<' '<<z<<' '<<curl[0]<<' '<<curl[1]<<' '<<curl[2]<<endl;
		  //std::cout<< "DIV B: "<<div<<endl;
		  //std::cout<< "CURL B: "<<curl[0]<<' '<<curl[1]<<' '<<curl[2]<<endl;
		  EXPECT_NEAR(curl[0], 0, 0.15);
		  EXPECT_NEAR(curl[1], 0, 0.15);
		  EXPECT_NEAR(curl[2], 0, 0.15);

	        }
    //fout.close();
    cout<<"bending radius: "<<myMagnet->Bend2D::designRadius_m<<endl;
    cout<<"field amplitude: "<<myMagnet->fieldAmplitude_m<<endl;
    /**
    Vector_t B(0.0);
    double phi_new = Physics::pi / 20 ;
    R[0] = (myMagnet->Bend2D::designRadius_m ) * cos(phi_new);
    R[1] = 0.005;
    R[2] = (myMagnet->Bend2D::designRadius_m ) * sin(phi_new);
    myMagnet->apply(R, P, 0 , E, B);
    cout<<"Derivative vs stepSize:"<<endl;
    cout<<partialsDerivB(R, B,1.0, myMagnet)[1][2]<<' '<<1.0<<endl;
    cout<<partialsDerivB(R, B,1.0e-1, myMagnet)[1][2]<<' '<<1.0e-1<<endl;
    cout<<partialsDerivB(R, B,1.0e-2, myMagnet)[1][2]<<' '<<1.0e-2<<endl;
    cout<<partialsDerivB(R, B,1.0e-3, myMagnet)[1][2]<<' '<<1.0e-3<<endl;
    cout<<partialsDerivB(R, B,1.0e-4, myMagnet)[1][2]<<' '<<1.0e-4<<endl;
    cout<<partialsDerivB(R, B,1.e-5, myMagnet)[1][2]<<' '<<1.0e-5<<endl;
    cout<<partialsDerivB(R, B,1.0e-6, myMagnet)[1][2]<<' '<<1.0e-6<<endl;
    cout<<partialsDerivB(R, B,1.0e-15, myMagnet)[1][2]<<' '<<1.0e-7<<endl;
    */
}

TEST(Quad, Quadrupole)
{
    OpalTestUtilities::SilenceTest silencer;

    MultipoleRep* quad = new MultipoleRep();
    //the following are used to initialise myMagnet
    PartData* partData = new PartData();
    PartBunch* bunch = new PartBunch(partData);
    bunch->resetM(0.938);
    bunch->setdT(1.0e-10);//time step
    double startField = 0.0, endField = 5.0 ;
    quad->setElementLength(3.0);
    ElementBase::ApertureType type = ElementBase::ApertureType::RECTANGULAR;
    vector<double> aperture(2, 1);
    quad->setAperture(type, aperture);
    quad->initialise(bunch, startField, endField);
    quad->setNormalComponent(2, 10.0);
    quad->setSkewComponent(2, 0.0);
    delete partData;
    delete bunch;
    Vector_t R(0. ,0., 0.);
    Vector_t P(0.), E(0.);
    ofstream gout("quadrupole_1");
    for(double z = -3.; z <= 9; z += 0.2)
    for(double x = -2; x <= 2; x += 0.01)
      for(double y = -10.0; y <= 10.; y += 1.)
        {
	  Vector_t B(0.0);
	  R(2) = z;
	  R(1) = y;
	  R(0) = x;
	  quad->apply(R, P, 0., E, B);
	  gout<<z<<' '<<x<<' '<<B[0]<<' '<<B[1]<<' '<<B[2]<<endl;
	  //gout<<x<<' '<<y<<' '<<sqrt(pow(B[0], 2.) + pow(B[1], 2.) + pow(B[2], 2.))<<endl;
	}
    gout.close();
    cout<<"length: "<<quad->getElementLength()<<endl;

}