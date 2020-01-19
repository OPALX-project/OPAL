// #include "Algorithms/Vektor.h"
// #include "Algorithms/PBunchDefs.h"
// #include "Algorithms/PartBunch.h"
// #include "Algorithms/PartData.h"
#include "Particle/IntCIC.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Particle/ParticleAttrib.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/Centering.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"

#include "opal_test_utilities/SilenceTest.h"

#include "gtest/gtest.h"
#define GUARDCELLSIZE 1
#define DIM 3

typedef IntCIC                                                IntrplCIC_t;

typedef Vektor<double, DIM>                                   Vector_t;
typedef ParticleSpatialLayout<double, DIM>::ParticlePos_t     Ppos_t;
typedef ParticleSpatialLayout<double, DIM>::ParticleIndex_t   PID_t;
typedef ParticleAttrib<double>                                Pscalar_t;
typedef InterpolatorTraits<double, DIM, IntrplCIC_t>::Cache_t Pcache_t;

typedef UniformCartesian<DIM>                                 Mesh_t;
typedef ParticleSpatialLayout< double, DIM, Mesh_t>           PLayout_t;

typedef Cell                                                  Edge_t;
typedef Cell                                                  Cell_t;
typedef Vert                                                  Vert_t;
typedef CenteredFieldLayout<DIM, Mesh_t, Edge_t>              FieldLayout_Edge_t;
typedef CenteredFieldLayout<DIM, Mesh_t, Cell_t>              FieldLayout_Cell_t;
typedef CenteredFieldLayout<DIM, Mesh_t, Vert_t>              FieldLayout_Vert_t;

typedef Field<Vector_t, DIM>                                  VField_Edge_t;
typedef Field<Vector_t, DIM>                                  VField_Cell_t;
typedef Field<Vector_t, DIM>                                  VField_Vert_t;

typedef Field<double, DIM, Mesh_t, Cell_t>                    SField_Cell_t;
typedef Field<double, DIM, Mesh_t, Vert_t>                    SField_Vert_t;
typedef GuardCellSizes<DIM>                                   GCS_t;
typedef NDIndex<DIM>                                          NDIdx_t;

typedef BConds<Vector_t, DIM, Mesh_t, Edge_t>                 BConds_t;
typedef ZeroGuardsAndZeroFace<Vector_t, DIM, Mesh_t, Edge_t>  GuardCell_Edge_t;
typedef PeriodicFace<Vector_t, DIM, Mesh_t, Edge_t>           PeriodicFace_Edge_t;

NDIdx_t getSingleElement(unsigned int i,
                         unsigned int j,
                         unsigned int k) {
    return NDIdx_t(Index(i,i),
                   Index(j,j),
                   Index(k,k));
}

class PartBunch: public IpplParticleBase<PLayout_t> {
public:
    PartBunch():
        IpplParticleBase<PLayout_t>()
    { }

    PartBunch(PLayout_t *pl):
        IpplParticleBase<PLayout_t>(pl)
    { }

    ParticleAttrib< double >     Q;
    ParticleAttrib< Vector_t >   P;
};

TEST(FieldLayoutTest, DISABLED_EdgeTest_1) {
    OpalTestUtilities::SilenceTest silencer;

    const unsigned int Nx = 32, Ny = 32, NZ = 64;
    e_dim_tag decomp[] = {PARALLEL, PARALLEL, SERIAL};
    double spacing[] = {1.0, 1.0, 1.0};

    // BConds_t vbc;
    // for (int i = 0; i < 6; ++ i) {
    //     vbc[i] = new PeriodicFace_Edge_t(i);
    // }

    Mesh_t mesh(Index(Nx),
                Index(Ny),
                Index(NZ),
                spacing,
                Vector_t(0.0));


    FieldLayout_Edge_t FL_edge(mesh, decomp);

    VField_Edge_t VFEdge(mesh,
                         FL_edge,
                         GuardCellSizes<DIM>(GUARDCELLSIZE));
    VFEdge = Vector_t(0.0);

    PLayout_t pl(FL_edge);
    PartBunch bunch(&pl);

    // bunch.create(1);

    // bunch.R[0] = Vector_t(0.8, 0.9, 0.7);
    // bunch.P[0] = Vector_t(1.0);

    // bunch.P.scatter(*VFEdge, bunch.R, IntCIC());

    // NDIdx_t elem = getSingleElement(0, 0, 0);
    // Vector_t value = VFEdge->localElement(elem);

    // EXPECT_NEAR(value[0], 0.8, 0.01);

    // std::cout << "EdgeCentering.cpp: " << __LINE__ << "\t" << value << std::endl;
}

/*
    Mesh_t *mesh_cc = new Mesh_t(Index(_Nx + 1),
                                 Index(_Ny + 1),
                                 Index(_Nz + 1),
                                 spacing,
                                 Vector_t(0.0));

    FieldLayout_Vert_t *FL_vert = new FieldLayout_Vert_t(*_mesh, decomp);

    SField_Vert_t *SFVert = new SField_Vert_t(*_mesh,
                                              *_FL_vert,
                                              //_vbc_edge,
                                              GuardCellSizes<3>(GUARDCELLSIZE));
    *SFVert = 0.0;

*/