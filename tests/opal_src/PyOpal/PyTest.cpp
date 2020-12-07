#include <stdlib.h>

#include "OpalSourcePath.h"
#include "gtest/gtest.h"

// not sure how to handle the path here
TEST(PyTest, RunPythonTests) {
  std::string opalSourceDir = OPAL_SOURCE_DIR; // tests/OpalSourcePath.h
  std::string sysCall = "python "+opalSourceDir+"/tests/opal_src/PyOpal/runTests.py";
  ASSERT_TRUE(system(NULL)) << "System calls not supported - failing";
  int i = system(sysCall.c_str());
  EXPECT_EQ(i, 0) << "Failed python unit tests";
}

