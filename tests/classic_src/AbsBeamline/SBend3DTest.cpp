/*
 *  Copyright (c) 2012, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <exception>
#include <stdexcept>
#include <fstream>
#include <cstdio>
#include <initializer_list>

#include "gtest/gtest.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/LogicalError.h"
#include "Structure/BoundaryGeometry.h"
#include "AbsBeamline/SBend3D.h"
#include "Physics/Physics.h"

#include "opal_test_utilities/SilenceTest.h"

#include "Utilities/LogicalError.h"

class LoadFieldMap {
  public:
    // Load a field map
    // name - file name to which the field map will be written
    // polynomial_order - order of the polynomial fit to use
    // smoothing_order - order at which the polynomial is smoothed (i.e.
    //                   derivatives set equal
    // skip_line - skip a line in the input file, index from 0; ignored if -ve
    LoadFieldMap(std::string name, int polynomial_order, int smoothing_order,
                 int skip_line)
        : sbend3d_m(nullptr), fname_m(name) {
        writeFieldMap(skip_line);
        getFieldMap(polynomial_order, smoothing_order);
    }

    ~LoadFieldMap() {
        delete sbend3d_m;
        removeFieldMap();
    }

    void getFieldMap(int polynomial_order, int smoothing_order) {
        try {
            sbend3d_m = new SBend3D("name");
            sbend3d_m->setLengthUnits(10.); // cm
            sbend3d_m->setFieldUnits(1.e-4); // ?
            sbend3d_m->setLengthUnits(10.); // cm
            sbend3d_m->setPolynomialOrder(polynomial_order); // ?
            sbend3d_m->setSmoothingOrder(smoothing_order); // ?
            sbend3d_m->setFieldMapFileName(fname_m);
          } catch (std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            sbend3d_m = nullptr;
            throw std::runtime_error("Failed to get field map");
        }
    }


    void writeFieldMap(int skip_line) {
        std::ofstream out(fname_m);
        std::initializer_list<std::string> init = {
           "      422280      422280      422280           1\n",
           " 1 X [LENGU]\n",
           "  2 Y [LENGU]\n",
           "  3 Z [LENGU]\n",
           "  4 BX [FLUXU]\n",
           "  5 BY [FLUXU]\n",
           "  6 BZ [FLUXU]\n",
           "  0\n",
           "    194.014700000       0.00000000000       80.3635200000      0.682759323466E-07  -5.37524925776      0.282807068056E-07\n",
           "    210.000000000       0.00000000000       0.00000000000       0.00000000000       5038.98826666       0.00000000000    \n",
           "    194.014700000       0.00000000000      -80.3635200000      0.682759055339E-07  -5.37524929545     -0.282807663112E-07\n",
           "    194.014700000      0.250000000000       80.3635200000     -0.250876890041E-01  -5.37468265486     -0.105656867773E-01\n",
           "    210.000000000      0.250000000000       0.00000000000       81.6064485044       5043.86454663     -0.162646432266E-02\n",
           "    194.014700000      0.250000000000      -80.3635200000     -0.252106804794E-01  -5.37468268550      0.102687594900E-01\n",
           "    194.014700000      0.500000000000       80.3635200000     -0.501213811240E-01  -5.37048637738     -0.211087378232E-01\n",
           "    210.000000000      0.500000000000       0.00000000000       163.479461107       5051.16103059     -0.319620022674E-02\n",
           "    194.014700000      0.500000000000      -80.3635200000     -0.503671731133E-01  -5.37048640042      0.205153440729E-01\n",
           "    194.938580000       0.00000000000       80.7462000000     -0.109193207295E-06  -5.47695895319     -0.452292335858E-07\n",
           "    211.000000000       0.00000000000       0.00000000000       0.00000000000       5346.80707376       0.00000000000    \n",
           "    194.938580000       0.00000000000      -80.7462000000     -0.109194150249E-06  -5.47695923342      0.452297706086E-07\n",
           "    194.938580000      0.250000000000       80.7462000000     -0.228949123764E-01  -5.47615607268     -0.965935126817E-02\n",
           "    211.000000000      0.250000000000       0.00000000000       71.6525328925       5351.82535885     -0.156196844700E-02\n",
           "    194.938580000      0.250000000000      -80.7462000000     -0.230215180683E-01  -5.47615608020      0.935369738637E-02\n",
           "    194.938580000      0.500000000000       80.7462000000     -0.457451523318E-01  -5.47168044969     -0.193017688931E-01\n",
           "    211.000000000      0.500000000000       0.00000000000       143.398642339       5359.42974721     -0.306846139443E-02\n",
           "    194.938580000      0.500000000000      -80.7462000000     -0.459981937280E-01  -5.47168045716      0.186908719298E-01\n",
           "    195.862460000       0.00000000000       81.1288900000       0.00000000000      -5.56838923159       0.00000000000    \n",
           "    212.000000000       0.00000000000       0.00000000000       0.00000000000       5612.41675566       0.00000000000    \n",
           "    195.862460000       0.00000000000      -81.1288900000       0.00000000000      -5.56838923767       0.00000000000    \n",
           "    195.862460000      0.250000000000       81.1288900000     -0.209034621315E-01  -5.56743204962     -0.884197597858E-02\n",
           "    212.000000000      0.250000000000       0.00000000000       61.0257302103       5617.50584022     -0.152992556601E-02\n",
           "    195.862460000      0.250000000000      -81.1288900000     -0.210330844408E-01  -5.56743205033      0.852904004824E-02\n",
           "    195.862460000      0.500000000000       81.1288900000     -0.416276950059E-01  -5.56278428219     -0.176093837788E-01\n",
           "    212.000000000      0.500000000000       0.00000000000       121.985754491       5625.18212699     -0.300403876153E-02\n",
           "    195.862460000      0.500000000000      -81.1288900000     -0.418867765797E-01  -5.56278428318      0.169839055390E-01\n"
        };
        std::vector<std::string> data(init);
        for (size_t i = 0; i < init.size(); ++i) {
            if (int(i) == skip_line)
                continue;
            out << data[i];
        }

        if (!out.good())
            throw std::runtime_error("could not write field map 'SBend3D_map.field'");

        out.close();
    }

    void removeFieldMap() {
        if (std::remove(fname_m.c_str()) != 0)
            std::cout << "Failed to delete field map " << fname_m
                      << " - it's probably okay" << std::endl;
    }


    SBend3D* sbend3d_m;
    std::string fname_m;
};

TEST(SBend3DTest, SBend3DGeometryTest) {
    OpalTestUtilities::SilenceTest silencer;

    LoadFieldMap fieldLoader("field1", 1, 1, -1);
    SBend3D* field = fieldLoader.sbend3d_m;
    if (field == nullptr)
        return; // skip the test
    Vector_t B, E, centroid;
    double radius = 2110.;
    for (double phi = -Physics::pi/4.+1e-3; phi < Physics::pi/2.; phi += Physics::pi/20)
        for (double r = -9.9; r < 31.; r += 5.)
            for (double y = -6.1; y < 6.; y += 1.) {
                Vector_t B(0., 0., 0.);
                Vector_t pos(
                    (radius+r)*cos(phi)-radius,
                    y,
                    (radius+r)*sin(phi)
                );
                field->apply(pos, Vector_t(0.0), 0, E, B);
                if (r > -10. && r < 10. &&
                    phi > 0. && phi < Physics::pi/4. &&
                    y > -5. && y < 5.) {
                    EXPECT_GT(fabs(B(1)), 1.e-4) // 1 Gauss
                          << "Pol:  " << r+radius << " " << y << " " << phi
                          << " Pos: " << pos(0) << " " << pos(1) << " " << pos(2)
                          << "   B: " << B(0) << " " << B(1) << " " << B(2)
                          << std::endl;
                } else {
                    EXPECT_LT(fabs(B(1)), 1.)
                          << "Pol:  " << r+radius << " " << y << " " << phi
                          << " Pos: " << pos(0) << " " << pos(1) << " " << pos(2)
                          << "   B: " << B(0) << " " << B(1) << " " << B(2)
                          << std::endl;
                }
    }
}

void testField(double r, double y, double phi,
               double bx, double by, double bz,
               double tol) {
    LoadFieldMap fieldLoader("field2", 1, 1, -1);
    SBend3D* field = fieldLoader.sbend3d_m;
    if (field == nullptr) {
        std::cerr << "SKIPPING SBEND3D TEST - FAILED TO LOAD FIELD" << std::endl;
        return; // skip the test
    }
    Vector_t B, E, centroid, pos;
    double radius = 2110.;
    pos = Vector_t(
        (radius+r)*cos(phi)-radius,
        y,
        (radius+r)*sin(phi)
    );
    field->apply(pos, Vector_t(0.0), 0, E, B);
    // the field map is rotated through pi/8. (into start at z=0. coordinates)
    double sR = sin(Physics::pi/8.);
    double cR = cos(Physics::pi/8.);
    double bxR = bx*cR-bz*sR;
    double bzR = bz*cR+bx*sR;
    EXPECT_NEAR(B(0), bxR, tol) << "R_p: " << r << " " << y << " " << phi
                                << " B: " << bx << " " << by << " " << bz;
    EXPECT_NEAR(B(1), by, tol) << "R_p: " << r << " " << y << " " << phi
                               << " B: " << bx << " " << by << " " << bz;
    EXPECT_NEAR(B(2), bzR, tol) << "R_p: " << r << " " << y << " " << phi
                                << " B: " << bx << " " << by << " " << bz;
}

TEST(SBend3DTest, SBend3DFieldTest) {
    OpalTestUtilities::SilenceTest silencer;

    // 211.000000000       0.00000000000       0.00000000000       0.00000000000              0.00000000000
    testField(0., 0., Physics::pi/8., 0., 5346.80707376*1e-4, 0., 1e-6);
    // 211.000000000      0.250000000000       0.00000000000       71.6525328925       5351.82535885     -0.156196844700E-02
    testField(0., 2.5, Physics::pi/8.,
               71.65253*1e-4, 5351.8253*1e-4, -0.1561E-02*1e-4, 1e-7);
}

TEST(SBend3DTest, SBend3DBadFileTest) {
    OpalTestUtilities::SilenceTest silencer;

    LoadFieldMap fieldLoader3("field5", 1, 1, -1); // should work okay
    EXPECT_THROW(
                 LoadFieldMap fieldLoader3("field6", 1, 1, 3), // missing header line
                 LogicalError
                 );
    EXPECT_THROW(
                 LoadFieldMap fieldLoader3("field7", 1, 1, 8), // missing first body
                 LogicalError
                 );
    EXPECT_THROW(
                 LoadFieldMap fieldLoader3("field8", 1, 1, 34), // missing last body
                 LogicalError
                 );
}

TEST(SBend3DTest, SBend3DPolyPatchTest) {
    OpalTestUtilities::SilenceTest silencer;

    // make the poly order > 1; this puts SBend3D in PolynomialPatch mode;
    // check the field map loads correctly including e.g. dipole symmetry
    try {
        LoadFieldMap fieldLoader1("field3", 2, 2, -1);
        LoadFieldMap fieldLoader2("field3", 2, 2, -1);
        LoadFieldMap fieldLoader3("field4", 1, 1, -1);
        SBend3D* field1 = fieldLoader1.sbend3d_m;
        SBend3D* field2 = fieldLoader2.sbend3d_m;
        SBend3D* field3 = fieldLoader3.sbend3d_m;
        const double dphi = Physics::pi/8;
        double radius = 2110.;
        for (double y = -5.; y < 5.1; y += 2.5)
            for (double r = 2100.; r < 2120.1; r += 5.)
                for (double phi =-dphi/2.; phi < dphi*3.1; phi += dphi/2.) {
                    Vector_t B, E, centroid;
                    Vector_t pos(r*cos(phi)-radius, y, r*sin(phi));
                    field1->apply(pos, Vector_t(0.0), 0, E, B);
                    field2->apply(pos, Vector_t(0.0), 0, E, B);
                    field3->apply(pos, Vector_t(0.0), 0, E, B);
                }
    } catch (LogicalError& err) {
        std::cerr << err.what() << std::endl;
    } catch (GeneralClassicException& err) {
        std::cerr << err.what() << std::endl;
    } catch (std::exception& err) {
        std::cerr << err.what() << std::endl;
    }
}


TEST(SBend3DTest, GeometryTest2) {
    OpalTestUtilities::SilenceTest silencer;

    // Sucked geometry information from
    //     Classic/AbsBeamline/Ring.cpp::appendElement
    // Transform in OPAL-T coords
    LoadFieldMap fieldLoader("field9", 1, 1, -1);
    SBend3D* field = fieldLoader.sbend3d_m;
    std::cerr << " SBend3DTest::GeometryTest2 A" << std::endl;
    Euclid3D delta = field->getGeometry().getTotalTransform();
    Vector3D v = delta.getVector();
    Vector3D r = delta.getRotation().getAxis();
    std::cerr << " SBend3DTest::GeometryTest2 B" << std::endl;

    // Transform to cycl coordinates
    Euclid3D euclid(v(2), v(0), -v(1), r(2), r(0), -r(1));
    delta = euclid;
    std::cerr << " SBend3DTest::GeometryTest2 C" << std::endl;
    // Calculate change in position
    Vector_t deltaPos(delta.getVector()(0), delta.getVector()(1), 0);
    double endRot = delta.getRotation().getAxis()(2);
    Vector_t deltaNorm(cos(endRot), sin(endRot), 0.);
    std::cerr << " SBend3DTest::GeometryTest2 D" << std::endl;

    std::cerr << 24.*(1-cos(M_PI/12.)) << " " << 24.*sin(M_PI/12.) << " ** " << cos(M_PI/12.) << " " << sin(M_PI/12.) << " ** " << M_PI/12. << std::endl;
    std::cerr << deltaPos << " ** " << deltaNorm << " ** " << endRot << std::endl;
}
