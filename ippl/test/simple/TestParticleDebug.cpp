// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit http://www.acl.lanl.gov/POOMS for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// TestParticleDebug.cpp , Tim Williams 8/11/1998
// This tests the functions [e,s]pap() and the function setFormat() from
// Utility/FieldDebug.[h,cpp] . These are meant to be called from the debugger,
// but this function tests whether they work (for a couple of possible calls)
// within a program. It also includes specialized function definitions like the
// user of ParticleDebug must have in his own source code in order to be able
// to access callable functions from the debugger, as an example for users.
// This function also tests the setInform() function, to specify the Inform
// object used internally by ParticleDebug functions. Constructing an Inform
// object that writes into a file makes it easy to do the comparson with
// correct output.

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Particle/IpplParticleBase.h"
#include "Particle/PAssign.h"
#include "Utility/RNGLattice.h"
#include "Utility/ParticleDebug.h"
#include "Utility/FieldDebug.h"

#ifdef IPPL_USE_STANDARD_HEADERS
#include <iostream>
#include <fstream>
using namespace std;
#else
#include <iostream>
#include <fstream>
#endif

//-----------------------------------------------------------------------------
// Simple user Particles class definition
class Particles: public IpplParticleBase< ParticleSpatialLayout<double, 3> > {
public:
  //tjwdebug: add a scalar attribute:
  ParticleAttrib<double> sa;
  // Constructor:
  Particles(ParticleSpatialLayout<double,3>* psl) :
    IpplParticleBase<ParticleSpatialLayout<double, 3> >(psl) {
      //tjwdebug: add a scalar attribute:
      addAttribute(sa);
  }
  // Destructor.
  virtual ~Particles() {}
  // Overload the = operator; does the same thing as the copy constructor.
  Particles& operator=(const Particles& p) {
    R = p.R;
    update();
    return(*this);
  }
};
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// User-inserted prototypes to get debugger access (examples for user ref.):
// Scalar (double):
void  dpap(ParticleAttrib<double>& pattr, bool docomm) {
  pap(pattr, docomm);
}
void depap(ParticleAttrib<double>& pattr, int i, bool docomm) {
  epap(pattr, i, docomm);
}
void dspap(ParticleAttrib<double>& pattr, int base, int bnd, int stride,
	   bool docomm) {
  spap(pattr, base, bnd, stride, docomm);
}
// 3D Vektor (double):
void  dv3pap(ParticleAttrib<Vektor<double,3> >& pattr, bool docomm) {
  pap(pattr, docomm);
}
void dv3epap(ParticleAttrib<Vektor<double,3> >& pattr, int i, bool docomm) {
  epap(pattr, i, docomm);
}
void dv3spap(ParticleAttrib<Vektor<double,3> >& pattr, int base, int bnd,
	     int stride, bool docomm) {
  spap(pattr, base, bnd, stride, docomm);
}
//-----------------------------------------------------------------------------

// Forward declarations:
void hardCodedOutput(char* filename); // Prototype of function defined below.
bool thediff(char* filename1, char* filename2);


int main(int argc, char *argv[])
{
  Ippl ippl(argc,argv);
  Inform testmsg(argv[0]);
  bool passed = true; // Pass/fail test

  // FieldLayout used for ParticleSpatialLayout below:
  int nx = 4, ny = 4, nz = 4;
  Index I(nx); Index J(ny); Index K(nz);
  // Specify multipple vnodes (8) to make sure this works right:
  FieldLayout<3> layout(I,J,K,PARALLEL,PARALLEL,PARALLEL,8);

  // Create Particles object
  ParticleSpatialLayout<double,3>* pslayout =
    new ParticleSpatialLayout<double,3>(layout);
  Particles parts(pslayout);
  int np = 16;
  parts.globalCreate(np);

  double deltaPX = 4.0/np;
  for (int p = 0; p < parts.getLocalNum(); p++) {
    double positionComponent = 0.0 + deltaPX*parts.ID[p];;
    parts.R[p](0) = positionComponent;
    parts.R[p](1) = positionComponent;
    parts.R[p](2) = positionComponent;
  }
  parts.update();

  // Inform output objects for test output:
  Inform* fdip =
    new Inform(NULL,"text.test.TestParticleDebug",Inform::OVERWRITE,0);
  Inform& fdi = *fdip;

  // --------------------------------------------------------------------------
  // WITH COMMUNICATION
  // --------------------------------------------------------------------------
  setPtclDbgInform(fdi);

  // Scalar ParticleAttribute -------------------------------------------------
  setFormat(16,1);

  fdi << endl << "--------pap(parts.ID, true)-------" << endl;
  pap(parts.ID, true);

  fdi << endl << "--------epap(parts.ID, np/4-1, true)-------" << endl;
  epap(parts.ID, np/4-1, true);

  fdi << endl << "--------spap(parts.ID, 0, np/4-1, 2, true)-------" << endl;
  spap(parts.ID, 0, np/4-1, 2, true);

  // 3D Vector ParticleAttribute ----------------------------------------------
  setFormat(1,8);

  fdi << endl << "--------pap(parts.R, true)-------" << endl;
  pap(parts.R, true);

  fdi << endl << "--------epap(parts.R, np/4-1, true)-------" << endl;
  epap(parts.R, np/4-1, true);

  fdi << endl << "--------spap(parts.R, 0, np/4-1, 2, true)-------" << endl;
  spap(parts.R, 0, np/4-1, 2, true);

  // --------------------------------------------------------------------------
  // NOW TURN OFF COMMUNICATION
  // --------------------------------------------------------------------------
  fdi.setPrintNode(INFORM_ALL_NODES);

  // Scalar ParticleAttribute -------------------------------------------------
  setFormat(16,1);

  fdi << endl << "--------pap(parts.ID, false)-------" << endl;
  pap(parts.ID, false);

  fdi << endl << "--------epap(parts.ID, np/4-1, false)-------" << endl;
  epap(parts.ID, np/4-1, false);

  fdi << endl << "--------spap(parts.ID, 0, np/4-1, 2, false)-------" << endl;
  spap(parts.ID, 0, np/4-1, 2, false);

  // 3D Vector ParticleAttribute ----------------------------------------------
  setFormat(1,8);

  fdi << endl << "--------pap(parts.R, false)-------" << endl;
  pap(parts.R, false);

  fdi << endl << "--------epap(parts.R, np/4-1, false)-------" << endl;
  epap(parts.R, np/4-1, false);

  fdi << endl << "--------spap(parts.R, 0, np/4-1, 2, false)-------" << endl;
  spap(parts.R, 0, np/4-1, 2, false);

  // Write out "by hand" into another file what the previous field-printing
  // functions should have produced; this will be compared with what they
  // actually did produce:
  hardCodedOutput("text.correct.TestParticleDebug");

  // Compare the two files by mocking up the Unix "diff" command:
  delete fdip;
  passed =
    thediff("text.test.TestParticleDebug","text.correct.TestParticleDebug");

  testmsg << ( (passed) ? "PASSED" : "FAILED" ) << endl;
  return 0;
}

//-----------------------------------------------------------------------------
// Mock up the Unix "diff" utility to compare two files:
//-----------------------------------------------------------------------------
bool thediff(char* filename1, char* filename2)
{
  bool same = true;
  char ch1, ch2;
  ifstream file1(filename1);
  ifstream file2(filename2);
  if (!file1) {
    cout << "thediff(): Failed to open file " << filename1 << endl;
    return(false);
  }
  if (!file2) {
    cout << "thediff(): Failed to open file " << filename2 << endl;
    return(false);
  }
  while (file1.get(ch1)) {          // Read file 1 char-by-char until eof
    if (file2.get(ch2)) {           // Read equivalent char from file 2
      if (ch1 != ch2) same = false; // If they're different,files are different
    }
    else {
      same = false;                 // If file 2 ends before file 1, different
    }
  }
  return same;
}

//-----------------------------------------------------------------------------
void hardCodedOutput(char* filename)
{
  ofstream of(filename);
  of << endl
     << "--------pap(parts.ID, true)-------" << endl
     << "....PE = 0 GLOBAL ptcle index subrange (0 : 15 : 1)...." << endl
     << "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 " << endl
     << "" << endl
     << "" << endl
     << "--------epap(parts.ID, np/4-1, true)-------" << endl
     << "....PE = 0 GLOBAL ptcle index subrange (3 : 3 : 1)...." << endl
     << "3 " << endl
     << "" << endl
     << "" << endl
     << "--------spap(parts.ID, 0, np/4-1, 2, true)-------" << endl
     << "....PE = 0 GLOBAL ptcle index subrange (0 : 2 : 2)...." << endl
     << "0 2 " << endl
     << "" << endl
     << "" << endl
     << "--------pap(parts.R, true)-------" << endl
     << "....PE = 0 GLOBAL ptcle index subrange (0 : 15 : 1)...." << endl
     << "( 0 , 0 , 0 ) " << endl
     << "( 0.25 , 0.25 , 0.25 ) " << endl
     << "( 0.5 , 0.5 , 0.5 ) " << endl
     << "( 0.75 , 0.75 , 0.75 ) " << endl
     << "( 1 , 1 , 1 ) " << endl
     << "( 1.25 , 1.25 , 1.25 ) " << endl
     << "( 1.5 , 1.5 , 1.5 ) " << endl
     << "( 1.75 , 1.75 , 1.75 ) " << endl
     << "( 2 , 2 , 2 ) " << endl
     << "( 2.25 , 2.25 , 2.25 ) " << endl
     << "( 2.5 , 2.5 , 2.5 ) " << endl
     << "( 2.75 , 2.75 , 2.75 ) " << endl
     << "( 3 , 3 , 3 ) " << endl
     << "( 3.25 , 3.25 , 3.25 ) " << endl
     << "( 3.5 , 3.5 , 3.5 ) " << endl
     << "( 3.75 , 3.75 , 3.75 ) " << endl
     << "" << endl
     << "" << endl
     << "--------epap(parts.R, np/4-1, true)-------" << endl
     << "....PE = 0 GLOBAL ptcle index subrange (3 : 3 : 1)...." << endl
     << "( 0.75 , 0.75 , 0.75 ) " << endl
     << "" << endl
     << "" << endl
     << "--------spap(parts.R, 0, np/4-1, 2, true)-------" << endl
     << "....PE = 0 GLOBAL ptcle index subrange (0 : 2 : 2)...." << endl
     << "( 0 , 0 , 0 ) " << endl
     << "( 0.5 , 0.5 , 0.5 ) " << endl
     << "" << endl
     << "" << endl
     << "--------pap(parts.ID, false)-------" << endl
     << "....PE = 0 LOCAL ptcle index range (0 : 15 : 1)...." << endl
     << "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 " << endl
     << "" << endl
     << "" << endl
     << "--------epap(parts.ID, np/4-1, false)-------" << endl
     << "....PE = 0 LOCAL ptcle index range (3 : 3 : 1)...." << endl
     << "3 " << endl
     << "" << endl
     << "" << endl
     << "--------spap(parts.ID, 0, np/4-1, 2, false)-------" << endl
     << "....PE = 0 LOCAL ptcle index range (0 : 3 : 2)...." << endl
     << "0 2 " << endl
     << "" << endl
     << "" << endl
     << "--------pap(parts.R, false)-------" << endl
     << "....PE = 0 LOCAL ptcle index range (0 : 15 : 1)...." << endl
     << "( 0 , 0 , 0 ) " << endl
     << "( 0.25 , 0.25 , 0.25 ) " << endl
     << "( 0.5 , 0.5 , 0.5 ) " << endl
     << "( 0.75 , 0.75 , 0.75 ) " << endl
     << "( 1 , 1 , 1 ) " << endl
     << "( 1.25 , 1.25 , 1.25 ) " << endl
     << "( 1.5 , 1.5 , 1.5 ) " << endl
     << "( 1.75 , 1.75 , 1.75 ) " << endl
     << "( 2 , 2 , 2 ) " << endl
     << "( 2.25 , 2.25 , 2.25 ) " << endl
     << "( 2.5 , 2.5 , 2.5 ) " << endl
     << "( 2.75 , 2.75 , 2.75 ) " << endl
     << "( 3 , 3 , 3 ) " << endl
     << "( 3.25 , 3.25 , 3.25 ) " << endl
     << "( 3.5 , 3.5 , 3.5 ) " << endl
     << "( 3.75 , 3.75 , 3.75 ) " << endl
     << "" << endl
     << "" << endl
     << "--------epap(parts.R, np/4-1, false)-------" << endl
     << "....PE = 0 LOCAL ptcle index range (3 : 3 : 1)...." << endl
     << "( 0.75 , 0.75 , 0.75 ) " << endl
     << "" << endl
     << "" << endl
     << "--------spap(parts.R, 0, np/4-1, 2, false)-------" << endl
     << "....PE = 0 LOCAL ptcle index range (0 : 3 : 2)...." << endl
     << "( 0 , 0 , 0 ) " << endl
     << "( 0.5 , 0.5 , 0.5 ) " << endl
     << "" << endl;
     of.close();
     return;
}

/***************************************************************************
 * $RCSfile: TestParticleDebug.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: TestParticleDebug.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $
 ***************************************************************************/
