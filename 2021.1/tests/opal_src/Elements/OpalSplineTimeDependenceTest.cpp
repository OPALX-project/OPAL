#include <sstream>

#include "gtest/gtest.h"
#include "Elements/OpalSplineTimeDependence.h"

#include "opal_test_utilities/SilenceTest.h"

TEST(OpalSplineTimeDependenceTest, ConstructorTest) {
    OpalTestUtilities::SilenceTest silencer;

    OpalSplineTimeDependence dep;
    OpalSplineTimeDependence* dep_clone = dep.clone("new name");
    EXPECT_EQ(dep_clone->getOpalName(), "new name");
}

TEST(OpalSplineSplineDependenceTest, PrintTest) {
    OpalTestUtilities::SilenceTest silencer;

    OpalSplineTimeDependence dep;
    std::stringstream _string;
    dep.print(_string);
    EXPECT_EQ(_string.str(), "SPLINE_TIME_DEPENDENCE;\n");
}
