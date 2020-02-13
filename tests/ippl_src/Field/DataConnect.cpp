#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "DataSource/DataConnect.h"
#include "DataSource/DataConnectCreator.h"
#include "Field/Field.h"
#include "FieldLayout/FieldLayout.h"
#include "Index/Index.h"

#include <iostream>

constexpr double roundOffError = 1e-10;

TEST(Field, DataConnect2D)
{
    OpalTestUtilities::SilenceTest silencer;
    int sizeX      = 100, sizeY   = 100;
    int centerX    =  5, centerY = 5;
    int iterations =  10;

    const unsigned Dim=2;
    Index I(sizeX);
    Index J(sizeY);
    Index PrintI(4);
    Index PrintJ(4);
    FieldLayout<Dim> layout(I,J,PARALLEL,PARALLEL,2*Ippl::getNodes());
    Field<double,Dim>  A(layout,GuardCellSizes<Dim>(1));
    Field<double,Dim> A2(layout,GuardCellSizes<Dim>(1));
    Field<double,Dim>  B(layout);

    // set up a DiscField for A
    DiscField<Dim> df("Ack", "Ack.out.config", 2, "Field<double,2>");

    // make other connections, using default connect method
    DataConnect *Bcon = DataConnectCreator::create("B Connection");
    B.connect("B", Bcon);

    A = 0.0;
    A[centerX][centerY] = 1.0*iterations;
    double fact = 1.0/9.0;

    for(int iter = 0 ; iter < iterations ; iter++ ) {
        std::cout << "Iteration " << iter << ": Computing new A ..." << std::endl;
        assign(B[I][J],
               fact*(A[I+1][J+1] + A[I+1][J  ] + A[I+1][J-1] +
                     A[I  ][J+1] + A[I  ][J  ] + A[I  ][J-1] +
                     A[I-1][J+1] + A[I-1][J  ] + A[I-1][J-1])
               );
        assign(A,B);
        std::cout << "  iter = " << iter << ", sum = " << sum(A) << std::endl;
        std::cout << "  new A = " << A[PrintI][PrintJ] << std::endl;
        std::cout << "  Writing A to file ..." << std::endl;
        df.write(A,0);
        B = -2.0;
        df.write(B,1);
        std::cout << "  Updating connection ..." << std::endl;
        B.updateConnection();
        if ((iter+1) % 10 == 0)
            A.interact();
    }

    EXPECT_NEAR(sum(A), 9.60864, 1e-5);

    // now read the data back into a Field
    std::cout << "....................................................." << std::endl;
    DiscField<Dim> readf("Ack", "Ack.in.config");
    unsigned int nrecords = readf.get_NumRecords();
    unsigned int nfields  = readf.get_NumFields();
    NDIndex<Dim> nsize    = readf.get_Domain();
    std::cout << "Reading back data, records = " << nrecords << ", fields = ";
    std::cout << nfields << ", total domain = " << nsize << std::endl;

    if (nfields > 0 && nrecords > 0 && nsize.size() > 0) {
        // read all records and print them out
        FieldLayout<Dim> rdlayout(nsize);
        Field<double,Dim> rdfield(rdlayout);
        for (unsigned int r=0; r < nrecords; ++r) {
            readf.read(rdfield, 0, r);
            std::cout << "Read record " << r << ": " << rdfield[PrintI][PrintJ] << std::endl;
        }
    } else {
        // there was an error reading the file
        std::cout << "An error was encountered while trying to read the file." <<std::endl;
        EXPECT_TRUE(false);
    }

    //EXPECT_NEAR(sum(readf), 9.60864, 1e-5);

    delete Bcon;
}

TEST(Field, DataConnect3D)
{
    OpalTestUtilities::SilenceTest silencer;
    int vnodes     = 16;
    int sizeX      = 20, sizeY   = 20, sizeZ   = 20;
    int centerX    = 10, centerY = 10, centerZ = 10;
    int iterations = 10;

    const unsigned Dim=3;
    Index I(sizeX);
    Index J(sizeY);
    Index K(sizeZ);
    FieldLayout<Dim> layout(I,J,K,PARALLEL,PARALLEL,PARALLEL, vnodes);
    Field<double,Dim> A(layout,GuardCellSizes<Dim>(1));
    Field<double,Dim> B(layout);
    DataConnect *conn = DataConnectCreator::create("Ack.config");
    A.connect("Ack");

    A = 0.0;
    A[centerX]  [centerY]  [centerZ]   = 512.0*iterations;
    A[centerX+3][centerY+7][0]         = 512.0*iterations;
    A[centerX-4][centerY-5][centerZ+8] = 128.0*iterations;
    A[sizeX-1]  [centerY+6][centerZ]   = 512.0*iterations;
    A[centerX+2][centerY]  [centerZ-5] = 512.0*iterations;

    double fact = 1.0/27.0;

    for(int iter = 0 ; iter < iterations ; iter++ ) {
        std::cout << "Computing new A ..." << std::endl;
        assign(B[I][J][K], fact*(A[I  ][J  ][K+1] + A[I  ][J  ][K-1] +
                                 A[I  ][J+1][K  ] + A[I  ][J-1][K  ] +
                                 A[I+1][J  ][K  ] + A[I-1][J  ][K  ]));
        assign(B[I][J][K], B[I][J][K] + fact*(A[I+1][J+1][K  ] + A[I+1][J-1][K  ] +
                                              A[I][J][K] +
                                              A[I-1][J+1][K  ] + A[I-1][J-1][K  ]));
        assign(B[I][J][K], B[I][J][K] + fact*(
                                              A[I+1][J+1][K+1] + A[I+1][J  ][K+1] + A[I+1][J-1][K+1] +
                                              A[I  ][J+1][K+1] +                    A[I  ][J-1][K+1] +
                                              A[I-1][J+1][K+1] + A[I-1][J  ][K+1] + A[I-1][J-1][K+1]));
        assign(B[I][J][K], B[I][J][K] + fact*(
                                              A[I+1][J+1][K-1] + A[I+1][J  ][K-1] + A[I+1][J-1][K-1] +
                                              A[I  ][J+1][K-1] +                    A[I  ][J-1][K-1] +
                                              A[I-1][J+1][K-1] + A[I-1][J  ][K-1] + A[I-1][J-1][K-1]));
        assign(A,B);
        std::cout << "iter = " << iter << ", sum = " << sum(A) << std::endl;
        std::cout << "Updating connection ..." << std::endl;
        A.updateConnection();
        std::cout << "iter = " << iter << ", sum = " << sum(A) << std::endl;
        A.interact();
    }

    EXPECT_NEAR(sum(A), 13235.1, 1e-1);

    delete conn;
}
