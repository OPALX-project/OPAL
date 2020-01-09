#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "AppTypes/Vektor.h"

constexpr unsigned int Dim = 3;
constexpr double margin = 1e-7;

TEST(Vektor, dot)
{
    OpalTestUtilities::SilenceTest silencer;

    //Vektor<double,Dim> boxMin(-1.0,-1.0,-1.0);
    Vektor<double,Dim> boxMax( 1.0, 1.0, 1.0);

    Vektor<double,Dim> p1( .5, .5, .5);
    Vektor<double,Dim> p2( 1.5, 1.5, 1.5);

    EXPECT_TRUE ( dot(p1,p1) <=  dot(boxMax,boxMax));
    EXPECT_FALSE( dot(p2,p2) <=  dot(boxMax,boxMax));

    EXPECT_NEAR ( dot(p1, boxMax), 1.5,  margin);
    EXPECT_NEAR ( dot(p2, p2),     6.75, margin);
}