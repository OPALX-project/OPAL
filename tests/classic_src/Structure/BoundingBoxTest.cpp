#include "gtest/gtest.h"

#include "Structure/BoundingBox.h"

#include "opal_test_utilities/SilenceTest.h"

TEST(BoundingBoxTest, InsideTest) {
    OpalTestUtilities::SilenceTest silencer;

    BoundingBox boundingBox = BoundingBox::getBoundingBox({Vector_t(0.0), Vector_t(1.0)});
    EXPECT_TRUE(boundingBox.isInside(Vector_t(0.5)));
    EXPECT_FALSE(boundingBox.isInside(Vector_t(1.01, 0.5, 0.5)));
    EXPECT_FALSE(boundingBox.isInside(Vector_t(0.5, 1.01, 0.5)));
    EXPECT_FALSE(boundingBox.isInside(Vector_t(0.5, 0.5, 1.01)));
}

TEST(BoundingBoxTest, IntersectionTest) {
    OpalTestUtilities::SilenceTest silencer;

    BoundingBox boundingBox = BoundingBox::getBoundingBox({Vector_t(0.0), Vector_t(1.0)});
    boost::optional<Vector_t> intersectionPoint = boundingBox.getIntersectionPoint(Vector_t(1.1, 0.5, 0.5), Vector_t(-1, 0, 0));
    EXPECT_TRUE(intersectionPoint.operator bool());
    EXPECT_NEAR(euclidean_norm(Vector_t(1.0, 0.5, 0.5) - intersectionPoint.value()), 0.0, 1e-8);

    intersectionPoint = boundingBox.getIntersectionPoint(Vector_t(1.1, 0.5, 0.5), Vector_t(1, 0, 0));
    EXPECT_FALSE(intersectionPoint.operator bool());

    intersectionPoint = boundingBox.getIntersectionPoint(Vector_t(1.1, 1.1, 0.5), Vector_t(-0.1, -1.0, 0));
    EXPECT_TRUE(intersectionPoint.operator bool());
    EXPECT_NEAR(euclidean_norm(Vector_t(1.0, 0.1, 0.5) - intersectionPoint.value()), 0.0, 1e-8);
}