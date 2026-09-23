//
// Class PartBunch
//   Particle Bunch.
//   A representation of a particle bunch as a vector of particles.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPAL_PartBunch_HH
#define OPAL_PartBunch_HH

#include "Algorithms/PartBunchBase.h"

#include <cstddef>
#include <functional>

class PartBunch: public PartBunchBase<double, 3> {

public:
    typedef IpplParticleBase<Layout_t> pbase_t;
    enum { Dim = Dimension };

public:

    /// Default constructor.
    //  Construct empty bunch.
    explicit PartBunch(const PartData *ref);

    PartBunch() = delete;
    PartBunch(const PartBunch &) = delete;
    PartBunch &operator=(const PartBunch &) = delete;

    ~PartBunch();

//     pbase_t* clone();

    void initialize(FieldLayout_t *fLayout);

    void do_binaryRepart();

    double getRho(int x, int y, int z);

    /*

      Mesh and Field Layout related functions

    */

    // MATTHIAS CHECK
    const Mesh_t &getMesh() const;

//     void setMesh(Mesh_t* mesh);

    // MATTHIAS CHECK
    Mesh_t &getMesh();

//     void setFieldLayout(FieldLayout_t* fLayout);

    // MATTHIAS CHECK
    FieldLayout_t &getFieldLayout();

    void setBCAllPeriodic();
    void setBCAllOpen();

    void setBCForDCBeam();

    /*
      Compatibility function push_back

    */

    VectorPair_t getEExtrema();

    void computeSelfFields();

    /** /brief used for self fields with binned distribution */
    void computeSelfFields(int b);

    void computeSelfFields_cycl(double gamma);
    void computeSelfFields_cycl(int b);

    /// Per-particle predicate on the local index, for scatterMasked() and gatherMasked().
    typedef std::function<bool(size_t)> ParticleMask_t;

    /** \brief CIC scatter of Q onto rho_m, of the particles with Bin >= 0 and accept(i) true.
     *
     * Without a predicate (an empty accept) it is Q.scatter(rho_m, R, IntrplCIC_t()) with Q
     * of the particles flagged lost (Bin < 0) set to zero for the call, which visits every
     * particle. With one, the others are not read, so they may lie outside the mesh, and a
     * predicate that is true for every particle gives the same rho_m, bit for bit.
     */
    void scatterMasked(const ParticleMask_t& accept);

    /** \brief CIC gather of eg_m into Ef, at the particles with accept(i) true.
     *
     * Ef of the others is left as it is, for the caller to set; they may lie outside the
     * mesh. A predicate that is true for every particle gives the Ef of
     * Ef.gather(eg_m, R, IntrplCIC_t()), bit for bit.
     */
    void gatherMasked(const ParticleMask_t& accept);

    void resetInterpolationCache(bool clearCache = false);

    void swap(unsigned int i, unsigned int j);

    /// scalar potential
    Field_t rho_m;

    /// vector field on the grid
    VField_t  eg_m;

    Inform &print(Inform &os);

private:

    void updateDomainLength(Vektor<int, 3>& grid);

    void updateFields(const Vector_t& hr, const Vector_t& origin);

    /// resize mesh to geometry specified
    void resizeMesh();

    /// for defining the boundary conditions
    BConds<double, 3, Mesh_t, Center_t> bc_m;
    BConds<Vector_t, 3, Mesh_t, Center_t> vbc_m;


    bool interpolationCacheSet_m;

    ParticleAttrib<CacheDataCIC<double, 3U> > interpolationCache_m;

    /// scatterMasked(): the charges and positions of the particles that a predicate accepts
    ParticleAttrib<double> maskedQ_m;
    ParticleAttrib<Vector_t> maskedR_m;

    //FIXME
    ParticleLayout<double, 3> & getLayout() {
        return pbase_m->getLayout();
    }

    //FIXME
    const ParticleLayout<double, 3>& getLayout() const {
        return pbase_m->getLayout();
    }
};



inline
double PartBunch::getRho(int x, int y, int z) {
    return rho_m[x][y][z].get();
}

inline
const Mesh_t &PartBunch::getMesh() const {
    const Layout_t* layout = static_cast<const Layout_t*>(&getLayout());
    return layout->getLayout().getMesh();
}

inline
Mesh_t &PartBunch::getMesh() {
    Layout_t* layout = static_cast<Layout_t*>(&getLayout());
    return layout->getLayout().getMesh();
}

// inline

inline Inform &operator<<(Inform &os, PartBunch &p) {
    return p.print(os);
}

#endif // OPAL_PartBunch_HH
