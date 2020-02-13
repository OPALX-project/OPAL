#include "gtest/gtest.h"

#include "opal_test_utilities/SilenceTest.h"

#include "AppTypes/Vektor.h"
#include "AppTypes/Tenzor.h"
#include "AppTypes/SymTenzor.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Index/Index.h"

constexpr unsigned int Dim = 2;
constexpr double margin = 1e-7;

TEST(AppTypes, Operators)
{
    OpalTestUtilities::SilenceTest silencer;
    const int size = 20;
    Index I(size), J(size);
    FieldLayout<Dim> layout(I,J);
    Field<double,Dim>
        F1(layout), F2(layout), F3(layout), F4(layout), F5(layout);
    Field< Vektor<double,Dim>, Dim, UniformCartesian<Dim> >
        V1(layout), V2(layout), V3(layout);
    Field< Tenzor<double,Dim>, Dim, UniformCartesian<Dim> >
        T1(layout), T2(layout), T3(layout);
    Field< SymTenzor<double,Dim>, Dim, UniformCartesian<Dim> >
        S1(layout), S2(layout), S3(layout);

    Vektor<double,Dim> Vinit1(1.0,2.0);
    Vektor<double,Dim> Vinit2(10.0,20.0);
    Vektor<double,Dim> Vinit3(100.0,200.0);

    Tenzor<double,Dim> Tinit1(1.0,2.0,2.0,3.0);
    Tenzor<double,Dim> Tinit2(10.0,20.0,20.0,30.0);
    Tenzor<double,Dim> Tinit3(100.0,200.0,200.0,300.0);

    SymTenzor<double,Dim> Sinit1(1.0,2.0,3.0);
    SymTenzor<double,Dim> Sinit2(10.0,20.0,30.0);
    SymTenzor<double,Dim> Sinit3(100.0,200.0,300.0);

    F1 = 1.0;
    F2 = 2.0;
    F3 = 3.0;
    V1 = Vinit1 ;
    V2 = Vinit2 ;
    V3 = Vinit3 ;
    T1 = Tinit1 ;
    T2 = Tinit2 ;
    T3 = Tinit3 ;
    S1 = Sinit1 ;
    S2 = Sinit2 ;
    S3 = Sinit3 ;

    // test Vektor Field operations

    V3 = V2 + V1;
    std::cout << " V3 " << std::endl << V3 << std::endl;
    EXPECT_NEAR(min(V3)[0], 11, margin);
    EXPECT_NEAR(max(V3)[1], 22, margin);
    V3 = V2 - V1;
    std::cout << " V3 " << std::endl << V3 << std::endl;
    EXPECT_NEAR(min(V3)[0],  9, margin);
    EXPECT_NEAR(max(V3)[1], 18, margin);
    V3 = F1 * V1;
    std::cout << " V3 " << std::endl << V3 << std::endl;
    EXPECT_NEAR(min(V3)[0], 1, margin);
    EXPECT_NEAR(max(V3)[1], 2, margin);
    V3 = V1 * F1;
    std::cout << " V3 " << std::endl << V3 << std::endl;
    EXPECT_NEAR(min(V3)[0], 1, margin);
    EXPECT_NEAR(max(V3)[1], 2, margin);
    F3 = dot(V1,V2);
    std::cout << " F3 " << std::endl << F3 << std::endl;
    EXPECT_NEAR(min(F3), 50 ,margin);

    // test Tenzor Field operations

    T3 = T2 + T1;
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], 11, margin);
    EXPECT_NEAR(max(T3)[1], 22, margin);
    EXPECT_NEAR(max(T3)[2], 22, margin);
    EXPECT_NEAR(max(T3)[3], 33, margin);
    T3 = T2 - T1;
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0],  9, margin);
    EXPECT_NEAR(max(T3)[1], 18, margin);
    EXPECT_NEAR(max(T3)[2], 18, margin);
    EXPECT_NEAR(max(T3)[3], 27, margin);
    T3 = dot(T1,T2);
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], 50, margin);
    EXPECT_NEAR(max(T3)[1], 80, margin);
    EXPECT_NEAR(max(T3)[2], 80, margin);
    EXPECT_NEAR(max(T3)[3],130, margin);
    F3 = dotdot(T1,T2);
    std::cout << " F3 " << std::endl << F3 << std::endl;
    EXPECT_NEAR(min(F3), 180 ,margin);
    V3 = dot(V1,T2);
    std::cout << " V3 " << std::endl << V3 << std::endl;
    EXPECT_NEAR(min(V3)[0], 50, margin);
    EXPECT_NEAR(max(V3)[1], 80, margin);
    V3 = dot(T2,V1);
    std::cout << " V3 " << std::endl << V3 << std::endl;
    EXPECT_NEAR(min(V3)[0], 50, margin);
    EXPECT_NEAR(max(V3)[1], 80, margin);
    T3 = F1 * T2;
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], 10, margin);
    EXPECT_NEAR(max(T3)[1], 20, margin);
    EXPECT_NEAR(max(T3)[2], 20, margin);
    EXPECT_NEAR(max(T3)[3], 30, margin);
    T3 = T2 * F1;
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], 10, margin);
    EXPECT_NEAR(max(T3)[1], 20, margin);
    EXPECT_NEAR(max(T3)[2], 20, margin);
    EXPECT_NEAR(max(T3)[3], 30, margin);

    // test SymTenzor Field operations

    S3 = S1 + S2;
    std::cout << " S3 " << std::endl << S3 << std::endl;
    EXPECT_NEAR(min(S3)[0], 11, margin);
    EXPECT_NEAR(max(S3)[1], 22, margin);
    EXPECT_NEAR(max(S3)[2], 33, margin);
    S3 = S1 - S2;
    std::cout << " S3 " << std::endl << S3 << std::endl;
    EXPECT_NEAR(min(S3)[0], -9, margin);
    EXPECT_NEAR(max(S3)[1],-18, margin);
    EXPECT_NEAR(max(S3)[2],-27, margin);
    S3 = dot(S1,S2);
    std::cout << " S3 " << std::endl << S3 << std::endl;
    EXPECT_NEAR(min(S3)[0], 50, margin);
    EXPECT_NEAR(max(S3)[1], 80, margin);
    EXPECT_NEAR(max(S3)[2],130, margin);
    F3 = dotdot(S1,S2);
    std::cout << " F3 " << std::endl << F3 << std::endl;
    EXPECT_NEAR(min(F3), 180 ,margin);
    V3 = dot(V1,S2);
    std::cout << " V3 " << std::endl << V3 << std::endl;
    EXPECT_NEAR(min(V3)[0], 50, margin);
    EXPECT_NEAR(max(V3)[1], 80, margin);
    V3 = dot(S2,V1);
    std::cout << " V3 " << std::endl << V3 << std::endl;
    EXPECT_NEAR(min(V3)[0], 50, margin);
    EXPECT_NEAR(max(V3)[1], 80, margin);
    T3 = T1 + S2;
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], 11, margin);
    EXPECT_NEAR(max(T3)[1], 22, margin);
    EXPECT_NEAR(max(T3)[2], 22, margin);
    EXPECT_NEAR(max(T3)[3], 33, margin);
    T3 = T1 - S2;
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], -9, margin);
    EXPECT_NEAR(max(T3)[1],-18, margin);
    EXPECT_NEAR(max(T3)[2],-18, margin);
    EXPECT_NEAR(max(T3)[3],-27, margin);
    T3 = dot(T1,S2);
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], 50, margin);
    EXPECT_NEAR(max(T3)[1], 80, margin);
    EXPECT_NEAR(max(T3)[2], 80, margin);
    EXPECT_NEAR(max(T3)[3],130, margin);
    F3 = dotdot(T1,S2);
    std::cout << " F3 " << std::endl << F3 << std::endl;
    EXPECT_NEAR(min(F3), 180 ,margin);
    T3 = S1 + T2;
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], 11, margin);
    EXPECT_NEAR(max(T3)[1], 22, margin);
    EXPECT_NEAR(max(T3)[2], 22, margin);
    EXPECT_NEAR(max(T3)[3], 33, margin);
    T3 = S2 - T1;
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0],  9, margin);
    EXPECT_NEAR(max(T3)[1], 18, margin);
    EXPECT_NEAR(max(T3)[2], 18, margin);
    EXPECT_NEAR(max(T3)[3], 27, margin);
    T3 = dot(S1,T2);
    std::cout << " T3 " << std::endl << T3 << std::endl;
    EXPECT_NEAR(min(T3)[0], 50, margin);
    EXPECT_NEAR(max(T3)[1], 80, margin);
    EXPECT_NEAR(max(T3)[2], 80, margin);
    EXPECT_NEAR(max(T3)[3],130, margin);
    F3 = dotdot(S1,T2);
    std::cout << " F3 " << std::endl << F3 << std::endl;
    EXPECT_NEAR(min(F3), 180 ,margin);
    S3 = F1 * S2;
    std::cout << " S3 " << std::endl << S3 << std::endl;
    EXPECT_NEAR(min(S3)[0], 10, margin);
    EXPECT_NEAR(max(S3)[1], 20, margin);
    EXPECT_NEAR(max(S3)[2], 30, margin);
    S3 = S2 * F1;
    std::cout << " S3 " << std::endl << S3 << std::endl;
    EXPECT_NEAR(min(S3)[0], 10, margin);
    EXPECT_NEAR(max(S3)[1], 20, margin);
    EXPECT_NEAR(max(S3)[2], 30, margin);
}