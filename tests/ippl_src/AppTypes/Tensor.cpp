#include "gtest/gtest.h"

#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"
#include "AppTypes/Tenzor.h"
#include "AppTypes/SymTenzor.h"
#include "AppTypes/AntiSymTenzor.h"

constexpr double margin = 1e-7;

template <int Dim>
class testTensorUnaryOps
{
public:
    static void apply()
    {
        NDIndex<Dim> domain;
        int nCellsTot = 1;
        for (int d = 0; d < Dim; d++) { 
            domain[d] = Index(5);
            nCellsTot *= 5;
        }
        FieldLayout<Dim> layout(domain);
        typedef UniformCartesian<Dim> Mesh;
        typedef Tenzor<double,Dim> FT_t;
        typedef AntiSymTenzor<double,Dim> AT_t;
        typedef SymTenzor<double,Dim> ST_t;
        Field<FT_t,Dim,Mesh> tff(layout);
        Field<AT_t,Dim,Mesh> tfa(layout);
        Field<ST_t,Dim,Mesh> tfs(layout);

        // Assign values:
        Tenzor<double,Dim> tf, tfTranspose;
        AntiSymTenzor<double,Dim> ta, taTranspose;
        SymTenzor<double,Dim> ts, tsTranspose;
        double fullSymTrace = 0.0;
        for (int i = 0; i < Dim; i++) {
            for (int j = 0; j < Dim; j++) {
                tf(i,j) = (i+1)*(i+1) + (j+1)*(j+1) + (i+4)*(j+4) + i;
                if (i == j) fullSymTrace += tf(i,j);
                tfTranspose(j,i) = tf(i,j);
            }
        }
        ta = tf;
        ts = tf;
        tff = tf;
        tfa = ta;
        tfs = ts;
        for (int i = 0; i < Dim; i++) {
            for (int j = 0; j < Dim; j++) {
                taTranspose(j,i) = ta(i,j);
                tsTranspose(j,i) = ts(i,j);
            }
        }

        // Test determinant of Tenzor:
        PInsist(Dim<4, "[Sym]Tenzor det() function not implemented for Dim>3!");
        double detValue = sum(det(tff));
        //  testmsg << "detValue = " << detValue << endl;
        switch (Dim) {
        case 1:
            EXPECT_NEAR(detValue, 18*nCellsTot, margin);
            break;
        case 2:
            EXPECT_NEAR(detValue, -38*nCellsTot, margin);
            break;
        case 3:
            EXPECT_NEAR(detValue, -4*nCellsTot, margin);
            break;
        default:
            ERRORMSG("Attempting to call det() for tensor greater than 3D" << endl);
            break;
        }

        // Test determinant of AntiSymTenzor
        double detValueA = sum(det(tfa));
        //  testmsg << "detValueA = " << detValueA << endl;
        switch (Dim) {
        case 1:
            EXPECT_NEAR(detValueA, 0, margin);
            break;
        case 2:
            EXPECT_NEAR(detValueA, -ta(1,0)*ta(0,1)*nCellsTot, margin);
            break;
        case 3:
            EXPECT_NEAR(detValueA, 0, margin);
            break;
        default:
            ERRORMSG("Attempting to call det() for tensor greater than 3D" << endl);
            break;
        }

        // Test determinant of SymTenzor
        double detValueS = sum(det(tfs));
        switch (Dim) {
        case 1:
            EXPECT_NEAR(detValueS, 18*nCellsTot, margin);
            break;
        case 2:
            EXPECT_NEAR(detValueS, -38.25*nCellsTot, margin);
            break;
        case 3:
            EXPECT_NEAR(detValueS, -4*nCellsTot, margin);
            break;
        default:
            ERRORMSG("Attempting to call det() for tensor greater than 3D" << endl);
            break;
        }

        // Test trace of Tenzor:
        double traceValue;
        traceValue = sum(trace(tff));
        EXPECT_NEAR(traceValue, fullSymTrace*nCellsTot, margin);

        // Test trace of AntiSymTenzor:
        traceValue = sum(trace(tfa));
        EXPECT_NEAR(traceValue, 0, margin);

        // Test trace of SymTenzor:
        traceValue = sum(trace(tfs));
        EXPECT_NEAR(traceValue, fullSymTrace*nCellsTot, margin);

        // Test transpose of Tenzor:
        Tenzor<double,Dim> transposeValue;
        transposeValue = sum(transpose(tff));
        EXPECT_TRUE((transposeValue == tfTranspose*nCellsTot));

        // Test transpose of AntiSymTenzor:
        AntiSymTenzor<double,Dim> transposeValueA;
        transposeValueA = sum(transpose(tfa));
        EXPECT_TRUE((transposeValueA == taTranspose*nCellsTot));

        // Test transpose of SymTenzor:
        SymTenzor<double,Dim> transposeValueS;
        transposeValueS = sum(transpose(tfs));
        EXPECT_TRUE((transposeValueS == tsTranspose*nCellsTot));

        // Test cofactors of Tenzor:
        Tenzor<double,Dim> cofactorsValue;
        cofactorsValue = sum(cofactors(tff))/nCellsTot;
        // Check results by computing det using all possible Laplace expansions,
        // and comparing to directly-computed det value:
        Tenzor<double,Dim> sumValue;
        sumValue = sum(tff);
        double altDetValue;
        // Laplace expansions using rows:
        for (int i = 0; i < Dim; i++) {
            altDetValue = 0.0;
            for (int j = 0; j < Dim; j++) {
                altDetValue += sumValue(i,j)*cofactorsValue(i,j);
            }
            EXPECT_NEAR(altDetValue, detValue, margin);
        }
        // Laplace expansions using columns:
        for (int j = 0; j < Dim; j++) {
            altDetValue = 0.0;
            for (int i = 0; i < Dim; i++) {
                altDetValue += sumValue(i,j)*cofactorsValue(i,j);
            }
            EXPECT_NEAR(altDetValue, detValue, margin);
        }

        // Test cofactors of AntiSymTenzor:
        cofactorsValue = sum(cofactors(tfa))/nCellsTot;
        // Check results by computing det using all possible Laplace expansions,
        // and comparing to directly-computed det value:
        sumValue = sum(tfa);
        // Laplace expansions using rows:
        for (int i = 0; i < Dim; i++) {
            altDetValue = 0.0;
            for (int j = 0; j < Dim; j++) {
                altDetValue += sumValue(i,j)*cofactorsValue(i,j);
            }
            EXPECT_NEAR(altDetValue, detValueA, margin);
        }
        // Laplace expansions using columns:
        for (int j = 0; j < Dim; j++) {
            altDetValue = 0.0;
            for (int i = 0; i < Dim; i++) {
                altDetValue += sumValue(i,j)*cofactorsValue(i,j);
            }
            EXPECT_NEAR(altDetValue, detValueA, margin);
        }

        // Test cofactors of SymTenzor:
        cofactorsValue = sum(cofactors(tfs))/nCellsTot;
        // Check results by computing det using all possible Laplace expansions,
        // and comparing to directly-computed det value:
        sumValue = sum(tfs);
        // Laplace expansions using rows:
        for (int i = 0; i < Dim; i++) {
            altDetValue = 0.0;
            for (int j = 0; j < Dim; j++) {
                altDetValue += sumValue(i,j)*cofactorsValue(i,j);
            }
            EXPECT_NEAR(altDetValue, detValueS, margin);
        }
        // Laplace expansions using columns:
        for (int j = 0; j < Dim; j++) {
            altDetValue = 0.0;
            for (int i = 0; i < Dim; i++) {
                altDetValue += sumValue(i,j)*cofactorsValue(i,j);
            }
            EXPECT_NEAR(altDetValue, detValueS, margin);
        }
    }
};

TEST(Tensor, Tensor)
{
    testTensorUnaryOps<1>::apply();
    testTensorUnaryOps<2>::apply();
    testTensorUnaryOps<3>::apply();
}