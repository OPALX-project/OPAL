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

// TestParticleDebugParallel.cpp , Tim Williams 8/11/1998
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

// JCC: Note that the definitions of the declared functions
// hardCodedOutput and thediff are not present, so this test
// code will not compile!

// include files
#include "Utility/IpplInfo.h"
#include "Utility/Inform.h"
#include "Message/Communicate.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Particle/IpplParticleBase.h"
#include "Particle/PAssign.h"
#include "Utility/RNGLattice.h"
#include "Utility/ParticleDebug.h"
#include "Utility/FieldDebug.h"

#include <iostream>

//-----------------------------------------------------------------------------
// Simple user Particles class definition
class Particles: public IpplParticleBase< ParticleSpatialLayout<double, 3> > {
public:
  // Constructor:
  Particles(ParticleSpatialLayout<double,3>* psl) :
    IpplParticleBase<ParticleSpatialLayout<double, 3> >(psl) {
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

  // FieldLayout used for ParticleSpatialLayout below:
  int nx = 4, ny = 4, nz = 4;
  Index I(nx); Index J(ny); Index K(nz);
  // Specify multipple vnodes (8) to make sure this works right:
  FieldLayout<3> layout(I,J,K,PARALLEL,PARALLEL,PARALLEL,8);

  // Create Particles object
  ParticleSpatialLayout<double, 3> pslayout(layout);
  Particles parts(&pslayout);
  int np = 16;
  parts.globalCreate(np);

  // Initialize position values
  parts.R(0) = RNGLattice<double> (0.0, 4.0, np);
  parts.R(1) = RNGLattice<double> (0.0, 4.0, np);
  parts.R(2) = RNGLattice<double> (0.0, 4.0, np);
  parts.update();

  // Inform output objects for test output:
  Inform fdi(NULL,INFORM_ALL_NODES);
  setPtclDbgInform(fdi);
  int mype = IpplInfo::myNode();

  // Scalar ParticleAttribute -------------------------------------------------
  setFormat(16,1);

  if (mype == 0) fdi << endl << "--------pap(parts.ID, true)-------" << endl;
  IpplInfo::Comm->barrier();
  pap(parts.ID, true);

  if (mype == 0)
    fdi << endl << "--------epap(parts.ID, np/4-1, true)-------" << endl;
  IpplInfo::Comm->barrier();
  epap(parts.ID, np/4-1, true);

  if (mype == 0)
    fdi << endl << "--------spap(parts.ID, 0, np/4-1, 2, true)-------" << endl;
  IpplInfo::Comm->barrier();
  spap(parts.ID, 0, np/4-1, 2, true);

  // 3D Vector ParticleAttribute ----------------------------------------------
  setFormat(1,8);

  if (mype == 0) fdi << endl << "--------pap(parts.R, true)-------" << endl;
  IpplInfo::Comm->barrier();
  pap(parts.R, true);

  if (mype == 0)
    fdi << endl << "--------epap(parts.R, np/4-1, true)-------" << endl;
  IpplInfo::Comm->barrier();
  epap(parts.R, np/4-1, true);

  if (mype == 0)
    fdi << endl << "--------spap(parts.R, 0, np/4-1, 2, true)-------" << endl;
  IpplInfo::Comm->barrier();
  spap(parts.R, 0, np/4-1, 2, true);

  // --------------------------------------------------------------------------
  // NOW TURN OFF COMMUNICATION
  // --------------------------------------------------------------------------

  // Scalar ParticleAttribute -------------------------------------------------
  setFormat(16,1);

  if (mype == 0) fdi << endl << "--------pap(parts.ID, false)-------" << endl;
  IpplInfo::Comm->barrier();
  int npes = IpplInfo::getNodes();
  for (int pe=0; pe<npes; pe++) {
    if (pe == mype) {
      pap(parts.ID, false);
    }
    IpplInfo::Comm->barrier();
  }

  if (mype == 0)
    fdi << endl << "--------epap(parts.ID, np/4-1, false)-------" << endl;
  IpplInfo::Comm->barrier();
  for (int pe=0; pe<npes; pe++) {
    if (pe == mype) {
      epap(parts.ID, np/4-1, false);
    }
    IpplInfo::Comm->barrier();
  }
  if (mype == 0)
    fdi << endl << "--------spap(parts.ID, 0, np/4-1, 2, false)-------"
	<< endl;
  IpplInfo::Comm->barrier();
  for (int pe=0; pe<npes; pe++) {
    if (pe == mype) {
      spap(parts.ID, 0, np/4-1, 2, false);
    }
    IpplInfo::Comm->barrier();
  }

  // 3D Vector ParticleAttribute ----------------------------------------------
  setFormat(1,8);

  if (mype == 0) fdi << endl << "--------pap(parts.R, false)-------" << endl;
  IpplInfo::Comm->barrier();
  for (int pe=0; pe<npes; pe++) {
    if (pe == mype) {
      pap(parts.R, false);
    }
    IpplInfo::Comm->barrier();
  }

  if (mype == 0)
    fdi << endl << "--------epap(parts.R, np/4-1, false)-------" << endl;
  IpplInfo::Comm->barrier();
  for (int pe=0; pe<npes; pe++) {
    if (pe == mype) {
      epap(parts.R, np/4-1, false);
    }
    IpplInfo::Comm->barrier();
  }

  if (mype == 0)
    fdi << endl << "--------spap(parts.R, 0, np/4-1, 2, false)-------" << endl;
  IpplInfo::Comm->barrier();
  for (int pe=0; pe<npes; pe++) {
    if (pe == mype) {
      spap(parts.R, 0, np/4-1, 2, false);
    }
    IpplInfo::Comm->barrier();
  }

  return 0;
}

/***************************************************************************
 * $RCSfile: TestParticleDebugParallel.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:39 $
 * IPPL_VERSION_ID: $Id: TestParticleDebugParallel.cpp,v 1.1.1.1 2003/01/23 07:40:39 adelmann Exp $
 ***************************************************************************/
