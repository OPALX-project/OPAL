#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "AppTypes/Vektor.h"
#include "Field/Field.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Meshes/CartesianCentering.h"
#include "Meshes/UniformCartesian.h"

#include <fstream>

// forward declarations
namespace {
    void hardCodedOutput(std::string filename); // Prototype of function defined below.
    bool thediff(std::string filename1, std::string filename2);

    CenteringEnum zz[2] = {CELL, VERTEX};
}

// template definitions
CenteringEnum CCCEnums<2U,1U,0U>::allCell[2U*1U];
CenteringEnum CCCEnums<3U,1U,1U>::allFace[3U*1U];
CenteringEnum CCCEnums<3U,3U,0U>::allVertex[3U*3U];
CenteringEnum CCCEnums<3U,3U,0U>::vectorFace[3U*3U];

TEST(Meshes, CartesianCentering)
{
    // For writing file output to compare against hardcoded correct file output:
    Inform fdi(NULL,"text.test.TestCartesianCentering",Inform::OVERWRITE,0);

    const unsigned nx=4, ny=4, nz=4;
    Index I(nx), J(ny), K(nz);

    const unsigned ND3 = 3;
    typedef UniformCartesian<ND3> M3;
    M3 m3(I,J,K);

    const unsigned ND2 = 2;
    typedef UniformCartesian<ND2> M2;
    M2 m2(I,J);

    typedef CommonCartesianCenterings<ND2,1U>::allCell CA;
    CenteredFieldLayout<ND2,M2,CA> clA(m2);
    Field<double, ND2, M2, CA > A(clA);
    A.print_Centerings(fdi.getStream());

    typedef CartesianCentering<zz,ND2,1U> CB;
    CenteredFieldLayout<ND2,M2,CB> clB(m2);
    Field<double, ND2, M2, CB > B(clB);
    B.print_Centerings(fdi.getStream());

    typedef CommonCartesianCenterings<ND3,1U,1U>::allFace CC;
    CenteredFieldLayout<ND3,M3,CC> clC(m3);
    Field<double, ND3, M3, CC> C(clC);
    C.print_Centerings(fdi.getStream());

    typedef CommonCartesianCenterings<ND3,3U>::allVertex CD;
    CenteredFieldLayout<ND3,M3,CD> clD(m3);
    Field<Vektor<double,ND3>, ND3, M3, CD> D(clD);
    D.print_Centerings(fdi.getStream());

    typedef CommonCartesianCenterings<ND3,3U>::vectorFace CE;
    CenteredFieldLayout<ND3,M3,CE> clE(m3);
    Field<Vektor<double,ND3>, ND3, M3, CE> E(clE);
    E.print_Centerings(fdi.getStream());

    fdi << endl ; // Needed to flush output to file

    // Write out "by hand" into another file what the previous field-printing
    // functions should have produced; this will be compared with what they
    // actually did produce:
    hardCodedOutput("text.correct.TestCartesianCentering");

    // Compare the two files by mocking up the Unix "diff" command:
    bool passed = thediff("text.test.TestCartesianCentering",
                          "text.correct.TestCartesianCentering");

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
        of << "CartesianCentering: no specialized name (yet) for this case" << std::endl
           << "Dim = 2 ; NComponents = 1" << std::endl
           << "centering[dim=0][component=0] = CELL  " << std::endl
           << "centering[dim=1][component=0] = CELL  " << std::endl
           << "CartesianCentering: no specialized name (yet) for this case" << std::endl
           << "Dim = 2 ; NComponents = 1" << std::endl
           << "centering[dim=0][component=0] = CELL  " << std::endl
           << "centering[dim=1][component=0] = VERTEX" << std::endl
           << "CartesianCentering: no specialized name (yet) for this case" << std::endl
           << "Dim = 3 ; NComponents = 1" << std::endl
           << "centering[dim=0][component=0] = CELL  " << std::endl
           << "centering[dim=1][component=0] = CELL  " << std::endl
           << "centering[dim=2][component=0] = CELL  " << std::endl
           << "CartesianCentering: no specialized name (yet) for this case" << std::endl
           << "Dim = 3 ; NComponents = 3" << std::endl
           << "centering[dim=0][component=0] = CELL  " << std::endl
           << "centering[dim=0][component=1] = CELL  " << std::endl
           << "centering[dim=0][component=2] = CELL  " << std::endl
           << "centering[dim=1][component=0] = CELL  " << std::endl
           << "centering[dim=1][component=1] = CELL  " << std::endl
           << "centering[dim=1][component=2] = CELL  " << std::endl
           << "centering[dim=2][component=0] = CELL  " << std::endl
           << "centering[dim=2][component=1] = CELL  " << std::endl
           << "centering[dim=2][component=2] = CELL  " << std::endl
           << "CartesianCentering: no specialized name (yet) for this case" << std::endl
           << "Dim = 3 ; NComponents = 3" << std::endl
           << "centering[dim=0][component=0] = CELL  " << std::endl
           << "centering[dim=0][component=1] = CELL  " << std::endl
           << "centering[dim=0][component=2] = CELL  " << std::endl
           << "centering[dim=1][component=0] = CELL  " << std::endl
           << "centering[dim=1][component=1] = CELL  " << std::endl
           << "centering[dim=1][component=2] = CELL  " << std::endl
           << "centering[dim=2][component=0] = CELL  " << std::endl
           << "centering[dim=2][component=1] = CELL  " << std::endl
           << "centering[dim=2][component=2] = CELL  " << std::endl
           << std::endl;
    of.close();
    return;
}
}