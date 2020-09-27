#ifndef PBUNCHDEFS_H
#define PBUNCHDEFS_H

#include "Algorithms/Vektor.h"
#include "Particle/IntCIC.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Meshes/UniformCartesian.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Field/Field.h"

#ifdef ENABLE_AMR
    #include "Amr/AmrDefs.h"
    #include "Amr/BoxLibParticle.h"
    #include "Amr/BoxLibLayout.h"
#endif

typedef IntCIC  IntrplCIC_t;

typedef ParticleSpatialLayout<double, 3>::ParticlePos_t Ppos_t;
typedef ParticleSpatialLayout<double, 3>::ParticleIndex_t PID_t;

typedef UniformCartesian<3, double> Mesh_t;

typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;

typedef Cell Center_t;

typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;

typedef Field<double, 3, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, 3, Mesh_t, Center_t>     VField_t;

#ifdef ENABLE_AMR
    typedef amr::AmrField_t                      AmrField_t;
    typedef amr::AmrScalarFieldContainer_t       AmrScalarFieldContainer_t;
    typedef amr::AmrVectorFieldContainer_t       AmrVectorFieldContainer_t;
    typedef BoxLibLayout<double, 3>              AmrLayout_t;
    typedef BoxLibParticle<AmrLayout_t>          AmrParticle_t;
#endif

#endif