/*
 *  Copyright (c) 2014, Chris Rogers
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
#include <fstream>

#include "gtest/gtest.h"

#include "opal_src/Utilities/MockComponent.h"

#include "BasicActions/DumpFields.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include "opal_test_utilities/SilenceTest.h"

#include "boost/filesystem.hpp"

namespace test {

void test() {}

void setOneAttribute(DumpFields* dump, std::string name, double value) {
    Attributes::setReal(*dump->findAttribute(name), value);
}

void setAttributes(DumpFields* dump,
                   double x0, double dx, double nx,
                   double y0, double dy, double ny,
                   double z0, double dz, double nz,
                   std::string filename) {
    setOneAttribute(dump, "X_START", x0);
    setOneAttribute(dump, "DX", dx);
    setOneAttribute(dump, "X_STEPS", nx);
    setOneAttribute(dump, "Y_START", y0);
    setOneAttribute(dump, "DY", dy);
    setOneAttribute(dump, "Y_STEPS", ny);
    setOneAttribute(dump, "Z_START", z0);
    setOneAttribute(dump, "DZ", dz);
    setOneAttribute(dump, "Z_STEPS", nz);
    Attributes::setString(*dump->findAttribute("FILE_NAME"), filename);
}

TEST(DumpFieldsTest, ConstructorDestructor) {
    OpalTestUtilities::SilenceTest silencer;

    // neither in the set and grid is null
    DumpFields* dump1 = new DumpFields();
    delete dump1;
    // grid is not null and it is in the set
    DumpFields* dump2 = new DumpFields();
    setAttributes(dump2, 1., 1., 1.,   1., 1., 1.,   1., 1., 1., "/dev/null");
    dump2->execute();
    delete dump2;
}

void execute_throws(DumpFields* dump, std::string reason) {
    try {
        dump->execute();
        EXPECT_TRUE(false) << reason;
    } catch (OpalException& exc) {
        // pass;
    }
}

TEST(DumpFieldsTest, executeTest) {
    OpalTestUtilities::SilenceTest silencer;

    // dump the fields
    DumpFields dump1;
    execute_throws(&dump1, "should throw due to nsteps < 1");
    setAttributes(&dump1, 1., 1., 1.,   1., 1., 1.,   1., 1., 1., "/dev/null");
    dump1.execute();  // should be okay (normal)
    setAttributes(&dump1, -1., -1., 1.,   -1., -1., 1.,   -1., -1., 1., "/dev/null");
    dump1.execute();  // should be okay (-ve step is okay)
    setAttributes(&dump1, -1., -1., 0.,   -1., -1., 1.,   -1., -1., 1., "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps x < 1");
    setAttributes(&dump1, -1., -1., 1.,   -1., -1., 0.,   -1., -1., 1., "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps y < 1");
    setAttributes(&dump1, -1., -1., 1.,   -1., -1., 1.,   -1., -1., 0., "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps z < 1");
    setAttributes(&dump1, -1., -1., 1.,   -1., -1., 1.,   -1., -1., 1.5, "/dev/null");
    execute_throws(&dump1, "should throw due to nsteps not integer");
}

void clear_files(std::set<std::string> const& files) {

    std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();

    for (const std::string& fname : files) {
        boost::filesystem::remove(Util::combineFilePath({auxDirectory, fname}));
    }
}

TEST(DumpFieldsTest, writeFieldsTest) {
    OpalTestUtilities::SilenceTest silencer;

    std::string auxDirectory = OpalData::getInstance()->getAuxiliaryOutputDirectory();

    boost::filesystem::create_directory(auxDirectory);

    std::string fname1 = "test1";
    std::string fname2 = "test2";
    std::string fname3 = "test3";
    std::string fname4 = "test4";

    clear_files({fname1, fname2, fname3, fname4});
    DumpFields dump1;
    setAttributes(&dump1, 1., 1., 1.,   1., 1., 1.,   1., 1., 1., fname1);
    dump1.execute();
    DumpFields dump2;
    setAttributes(&dump2, 1., 1., 1.,   1., 1., 1.,   1., 1., 1., fname2);
    dump2.execute();
    DumpFields dump3;
    setAttributes(&dump3, 1., 1., 1.,   1., 1., 1.,   1., 1., 1., fname3);
    // note we don't execute dump3; so it should not be written
    DumpFields dump4;
    setAttributes(&dump4, 0.1, 0.1, 3.,   -0.1, 0.2, 2.,   0.2, 0.3, 2., fname4);
    dump4.execute();
    MockComponent comp;
    try {
        DumpFields::writeFields(&comp);
    } catch (OpalException& exc) {
        EXPECT_TRUE(false) << "Threw OpalException on writefields: " << exc.what() << std::endl;;
    }
    std::ifstream fin1(Util::combineFilePath({auxDirectory, fname1}));
    EXPECT_TRUE(fin1.good());
    std::ifstream fin2(Util::combineFilePath({auxDirectory, fname2}));
    EXPECT_TRUE(fin2.good());
    std::ifstream fin3(Util::combineFilePath({auxDirectory, fname3}));
    EXPECT_FALSE(fin3.good());  // does not exist
    std::ifstream fin4(Util::combineFilePath({auxDirectory, fname4}));
    EXPECT_TRUE(fin4.good());
    int n_lines;
    fin4 >> n_lines;
    EXPECT_EQ(n_lines, 12);
    std::string test_line;
    for (size_t i = 0; i < 8; ++i) {
        std::getline(fin4, test_line);
        // std::cerr << test_line << std::endl;
    }
    std::vector<double> line(6, 0.);
    for (size_t line_index = 0; line_index < 12; ++line_index) {
        for (size_t i = 0; i < 6; ++i) {
            fin4 >> line[i];
        }
        if (line_index == 0) {
            ASSERT_EQ(line[0], 0.1);
            ASSERT_EQ(line[1], -0.1);
            ASSERT_EQ(line[2], 0.2);
        }
        if (line[1] < 0.) {
            ASSERT_EQ(line[3], line[0]);
            ASSERT_EQ(line[4], line[1]);
            ASSERT_EQ(line[5], line[2]);
        } else {
            ASSERT_EQ(line[3], 0.);
            ASSERT_EQ(line[4], 0.);
            ASSERT_EQ(line[5], 0.);
        }
    }
    clear_files({fname1, fname2, fname3, fname4});
}

}
