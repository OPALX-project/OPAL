#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "FieldLayout/FieldLayout.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Particle/IpplParticleBase.h"
#include "Particle/PAssign.h"
#include "Utility/ParticleDebug.h"
#include "Utility/FieldDebug.h"

#include <fstream>
#include <iostream>

// TestParticleDebug , Tim Williams 8/11/1998
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

namespace {
    void hardCodedOutput(std::string filename); // Prototype of function defined below.
    bool thediff(std::string filename1, std::string filename2);
}

TEST(Particle, ParticleDebug)
{
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
    for (unsigned int p = 0; p < parts.getLocalNum(); p++) {
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
    bool passed =
        thediff("text.test.TestParticleDebug","text.correct.TestParticleDebug");
    EXPECT_TRUE(passed);
}

namespace {
    //-----------------------------------------------------------------------------
    // Mock up the Unix "diff" utility to compare two files:
    //-----------------------------------------------------------------------------
    bool thediff(std::string filename1, std::string filename2)
    {
        bool same = true;
        char ch1, ch2;
        std::ifstream file1(filename1);
        std::ifstream file2(filename2);
        if (!file1) {
            std::cout << "thediff(): Failed to open file " << filename1 << std::endl;
            return(false);
        }
        if (!file2) {
            std::cout << "thediff(): Failed to open file " << filename2 << std::endl;
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
    void hardCodedOutput(std::string filename)
    {
        std::ofstream of(filename);
        of << std::endl
           << "--------pap(parts.ID, true)-------" << std::endl
           << "....PE = 0 GLOBAL ptcle index subrange (0 : 15 : 1)...." << std::endl
           << "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------epap(parts.ID, np/4-1, true)-------" << std::endl
           << "....PE = 0 GLOBAL ptcle index subrange (3 : 3 : 1)...." << std::endl
           << "3 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------spap(parts.ID, 0, np/4-1, 2, true)-------" << std::endl
           << "....PE = 0 GLOBAL ptcle index subrange (0 : 2 : 2)...." << std::endl
           << "0 2 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------pap(parts.R, true)-------" << std::endl
           << "....PE = 0 GLOBAL ptcle index subrange (0 : 15 : 1)...." << std::endl
           << "( 0 , 0 , 0 ) " << std::endl
           << "( 0.25 , 0.25 , 0.25 ) " << std::endl
           << "( 0.5 , 0.5 , 0.5 ) " << std::endl
           << "( 0.75 , 0.75 , 0.75 ) " << std::endl
           << "( 1 , 1 , 1 ) " << std::endl
           << "( 1.25 , 1.25 , 1.25 ) " << std::endl
           << "( 1.5 , 1.5 , 1.5 ) " << std::endl
           << "( 1.75 , 1.75 , 1.75 ) " << std::endl
           << "( 2 , 2 , 2 ) " << std::endl
           << "( 2.25 , 2.25 , 2.25 ) " << std::endl
           << "( 2.5 , 2.5 , 2.5 ) " << std::endl
           << "( 2.75 , 2.75 , 2.75 ) " << std::endl
           << "( 3 , 3 , 3 ) " << std::endl
           << "( 3.25 , 3.25 , 3.25 ) " << std::endl
           << "( 3.5 , 3.5 , 3.5 ) " << std::endl
           << "( 3.75 , 3.75 , 3.75 ) " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------epap(parts.R, np/4-1, true)-------" << std::endl
           << "....PE = 0 GLOBAL ptcle index subrange (3 : 3 : 1)...." << std::endl
           << "( 0.75 , 0.75 , 0.75 ) " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------spap(parts.R, 0, np/4-1, 2, true)-------" << std::endl
           << "....PE = 0 GLOBAL ptcle index subrange (0 : 2 : 2)...." << std::endl
           << "( 0 , 0 , 0 ) " << std::endl
           << "( 0.5 , 0.5 , 0.5 ) " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------pap(parts.ID, false)-------" << std::endl
           << "....PE = 0 LOCAL ptcle index range (0 : 15 : 1)...." << std::endl
           << "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------epap(parts.ID, np/4-1, false)-------" << std::endl
           << "....PE = 0 LOCAL ptcle index range (3 : 3 : 1)...." << std::endl
           << "3 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------spap(parts.ID, 0, np/4-1, 2, false)-------" << std::endl
           << "....PE = 0 LOCAL ptcle index range (0 : 3 : 2)...." << std::endl
           << "0 2 " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------pap(parts.R, false)-------" << std::endl
           << "....PE = 0 LOCAL ptcle index range (0 : 15 : 1)...." << std::endl
           << "( 0 , 0 , 0 ) " << std::endl
           << "( 0.25 , 0.25 , 0.25 ) " << std::endl
           << "( 0.5 , 0.5 , 0.5 ) " << std::endl
           << "( 0.75 , 0.75 , 0.75 ) " << std::endl
           << "( 1 , 1 , 1 ) " << std::endl
           << "( 1.25 , 1.25 , 1.25 ) " << std::endl
           << "( 1.5 , 1.5 , 1.5 ) " << std::endl
           << "( 1.75 , 1.75 , 1.75 ) " << std::endl
           << "( 2 , 2 , 2 ) " << std::endl
           << "( 2.25 , 2.25 , 2.25 ) " << std::endl
           << "( 2.5 , 2.5 , 2.5 ) " << std::endl
           << "( 2.75 , 2.75 , 2.75 ) " << std::endl
           << "( 3 , 3 , 3 ) " << std::endl
           << "( 3.25 , 3.25 , 3.25 ) " << std::endl
           << "( 3.5 , 3.5 , 3.5 ) " << std::endl
           << "( 3.75 , 3.75 , 3.75 ) " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------epap(parts.R, np/4-1, false)-------" << std::endl
           << "....PE = 0 LOCAL ptcle index range (3 : 3 : 1)...." << std::endl
           << "( 0.75 , 0.75 , 0.75 ) " << std::endl
           << "" << std::endl
           << "" << std::endl
           << "--------spap(parts.R, 0, np/4-1, 2, false)-------" << std::endl
           << "....PE = 0 LOCAL ptcle index range (0 : 3 : 2)...." << std::endl
           << "( 0 , 0 , 0 ) " << std::endl
           << "( 0.5 , 0.5 , 0.5 ) " << std::endl
           << "" << std::endl;
        of.close();
        return;
    }
}