#include "gtest/gtest.h"

#include "Index/Index.h"
#include "FieldLayout/FieldLayout.h"
#include "Field/Field.h"
#include "Meshes/UniformCartesian.h"
#include "AppTypes/Vektor.h"

#include <cmath>

const unsigned Dim = 2;
int size = 8;
constexpr double margin = 1e-7;
typedef Vektor<double,Dim> Vek;

    Vek vabs(Vek v1)
    {
        Vek v2(v1);
        for (unsigned d=0; d<Dim; ++d)
            if ( v2[d]<0 )
                v2[d] = -v2[d];
        return v2;
    }

    Vek vmax(Vek v1, Vek v2)
    {
        Vek v0;
        for (unsigned d=0; d<Dim; ++d)
            v0[d] = std::max(v1[d],v2[d]);
        return v0;
    }

    double cutoff;
    Vek vselect1(Vek v1, Vek v2)
    {
        Vek v0 = vabs(v1 - v2);
        for (unsigned d=0; d<Dim; ++d)
            if ( v0[d] > cutoff )
                v0[d] = v2[d];
            else
                v0[d] = v1[d];
        return v0;
    }

    Vek vselect2(Vek diff, Vek cutoff)
    {
        Vek v0;
        for (unsigned d=0; d<Dim; ++d)
            if ( cutoff[d] >= fabs(diff[d])  )
                v0[d] = diff[d];
            else
                v0[d] = 0.0;
        return v0;
    }

    UNARY_FUNCTION(Vek,vabs,Vek)
    BINARY_FUNCTION(Vek,vmax,Vek, Vek)
    BINARY_FUNCTION(Vek,vselect1,Vek,Vek)
    BINARY_FUNCTION(Vek,vselect2,Vek,Vek)


TEST(Vektor, tz)
{
    Index I(size), J(size);
    int i,j;
    FieldLayout<Dim> layout(I,J);
    typedef Field< Vek ,Dim,UniformCartesian<Dim> > FV;
    typedef Field< double ,Dim,UniformCartesian<Dim> > Fd;
    FV V1(layout),V2(layout),V3(layout);
    Fd d1(layout),d2(layout),d3(layout);
    double err;

    for (i=0; i<size; ++i)
        for (j=0; j<size; ++j) {
            V1[i][j] = Vek( double(i), double(j) ) ;
            V2[i][j] = Vek( double(j), double(i) ) ;
        }

    V3 = 0.0;
    assign( V3, vabs(V2-V1)+V3+vabs(V3) );

    err = 0;
    for (i=0; i<size; ++i)
        for (j=0; j<size; ++j)
            {
                Vek v = V3[i][j].get();
                v -= Vek( abs(j-i), abs(i-j) );
                err += fabs(v[0]);
                err += fabs(v[1]);
            }
    EXPECT_NEAR(err,0,margin);

    assign( V3, vmax(V2,V1) );
    err = 0;
    for (i=0; i<size; ++i)
        for (j=0; j<size; ++j)
            {
                Vek v = V3[i][j].get();
                double m = (i>j) ? i : j;
                v -= Vek( m,m );
                err += fabs(v[0]);
                err += fabs(v[1]);
            }
    EXPECT_NEAR(err,0,margin);

    cutoff = 1.0;
    assign( V3, vselect1(V1,V2) );
    err = 0;
    for (i=0; i<size; ++i)
        for (j=0; j<size; ++j)
            {
                Vek v = V3[i][j].get();
                v -= vselect1( Vek(i,j), Vek(j,i) );
                err += fabs(v[0]);
                err += fabs(v[1]);
            }
    EXPECT_NEAR(err,0,margin);

    V3 = 1.0;
    V2 += vselect2(V1-V2,V3) ;
    err = 0;
    for (i=0; i<size; ++i)
        for (j=0; j<size; ++j)
            {
                Vek v = V2[i][j].get();
                v -= vselect1( Vek(i,j), Vek(j,i) );
                err += fabs(v[0]);
                err += fabs(v[1]);
            }
    EXPECT_NEAR(err,0,margin);
}