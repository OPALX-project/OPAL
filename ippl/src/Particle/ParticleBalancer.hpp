// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Particle/ParticleBalancer.h"
#include "Particle/IpplParticleBase.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Particle/ParticleUniformLayout.h"
#include "Particle/ParticleAttrib.h"
#include "Particle/IntNGP.h"
#include "Region/RegionLayout.h"
#include "Index/NDIndex.h"
#include "FieldLayout/FieldLayout.h"
#include "FieldLayout/BinaryBalancer.h"
#include "Utility/IpplInfo.h"



/////////////////////////////////////////////////////////////////////////////
// number density scatter as scatter(f, pp, intop), but only of the particles
// inside the given region.  A particle outside it is not counted, and its
// position is never looked up on the Field: beyond the guard cells, the lookup
// would abort.
template <class FT, unsigned Dim, class M, class C, class PT, class IntOp>
void
scatterInsideDomain(Field<FT,Dim,M,C>& f,
                    const ParticleAttrib< Vektor<PT,Dim> >& pp,
                    const NDRegion<PT,Dim>& domain, const IntOp& /*intop*/) {

  // make sure field is uncompressed and guard cells are zeroed
  f.Uncompress();
  FT zero = 0;
  f.setGuardCells(zero);

  const M& mesh = f.get_mesh();
  // iterate through the particles and call scatter operation for those
  // inside the half-open domain [min, max)
  typename ParticleAttrib< Vektor<PT,Dim> >::const_iterator ppiter;
  size_t i = 0;
  for (ppiter = pp.cbegin(); i < pp.size(); ++i, ++ppiter) {
    const Vektor<PT,Dim>& pos = *ppiter;
    bool inside = true;
    for (unsigned int d = 0; d < Dim; ++d)
      inside = inside && (pos[d] >= domain[d].min() && pos[d] < domain[d].max());
    if (inside)
      IntOp::scatter(FT(1), f, pos, mesh);
  }

  // accumulate values in guard cells
  f.accumGuardCells();

  INCIPPLSTAT(incParticleScatters);
}


/////////////////////////////////////////////////////////////////////////////
// calculate a new RegionLayout for a given ParticleBase, and distribute the
// new RegionLayout to all the nodes.  This uses a Field BinaryBalancer.
// With ParticleSpatialLayout::setOutsideToNearest(true), the particles
// outside the domain of the RegionLayout do not count in the balance.
template < class T, unsigned Dim, class Mesh, class CachingPolicy>
bool
BinaryRepartition(IpplParticleBase<ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy> >& PB, double offset) {



  static IntNGP interp; // to scatter particle density

  //Inform dbgmsg("Particle BinaryRepartition", INFORM_ALL_NODES);
  //dbgmsg << "Performing particle load balancing, for ";
  //dbgmsg << PB.getTotalNum() << " particles ..." << endl;

  // get the internal FieldLayout from the Particle object's internal
  // RegionLayout.  From this, we make a new Field (we do not need a
  RegionLayout<T,Dim,Mesh>& RL = PB.getLayout().getLayout();
  if ( ! RL.initialized()) {
    ERRORMSG("Cannot repartition particles: uninitialized layout." << endl);
    return false;
  }
  FieldLayout<Dim>& FL = RL.getFieldLayout();
  Mesh& mesh = RL.getMesh();

  // NDIndex which describes the entire domain ... if a particle is
  // outside this region, we are in trouble!
  const NDIndex<Dim>& TotalDomain = FL.getDomain();

  // for all the particles, do the following:
  //   1. get the position, and invert to the 'FieldLayout' index space
  //   2. increment the field at the position near this index position
  NDIndex<Dim> indx;

  // By default, we do the number density computation and repartition of
  // index space on a cell-centered Field.  If FieldLayout is vertex-centered,
  // we'll need to make some adjustments here.
  bool CenterOffset[Dim];
  int CenteringTotal = 0;
  unsigned int d;
  for (d=0; d<Dim; ++d) {
    CenterOffset[d] = (TotalDomain[d].length() < mesh.gridSizes[d]);
    CenteringTotal += CenterOffset[d];
  }


  // particles outside the domain are legal with this setting
  const bool outsideToNearest = PB.getLayout().getOutsideToNearest();

  if (CenteringTotal == Dim) { // allCell centering
    Field<double,Dim,Mesh,Cell> BF(mesh,FL,GuardCellSizes<Dim>(1));

    // Now do a number density scatter on this Field
    // Afterwards, the Field will be deleted, and will checkout of the
    // FieldLayout.  This is desired so that when we repartition the
    // FieldLayout, we do not waste time redistributing the Field's data.
    BF = offset;
    if (outsideToNearest)
      scatterInsideDomain(BF,PB.R,RL.getDomain(),interp);
    else
      scatter(BF,PB.R,interp);

    // calculate a new repartitioning of the field, and use this to repartition
    // the FieldLayout used inside the Particle object
    try
	{
		indx = CalcBinaryRepartition(FL, BF);
	}
	catch(BinaryRepartitionFailed bf)
	{
		return false;
	}
  }
  else if (CenteringTotal == 0) { // allVert centering
    Field<double,Dim,Mesh,Vert> BF(mesh,FL,GuardCellSizes<Dim>(1));

    // Now do a number density scatter on this Field
    // Afterwards, the Field will be deleted, and will checkout of the
    // FieldLayout.  This is desired so that when we repartition the
    // FieldLayout, we do not waste time redistributing the Field's data.
    BF = offset;
    if (outsideToNearest)
      scatterInsideDomain(BF,PB.R,RL.getDomain(),interp);
    else
      scatter(BF,PB.R,interp);

    // calculate a new repartitioning of the field, and use this to repartition
    // the FieldLayout used inside the Particle object
	try
	{
		indx = CalcBinaryRepartition(FL, BF);
	}
	catch(BinaryRepartitionFailed bf)
	{
		return false;
	}
  }
  else {
    ERRORMSG("Not implemented for face- and edge-centered Fields!!" << endl);
    Ippl::abort();
  }

  // now, we can repartition the FieldLayout within the RegionLayout
  RL.RepartitionLayout(indx);
  PB.update();
  return true;
}


// the same, but taking a uniform layout (this will not actually do anything)
template<class T, unsigned Dim>
bool
BinaryRepartition(IpplParticleBase<ParticleUniformLayout<T,Dim> >& /*PB*/, double /*offset*/) {
  // for a uniform layout, this repartition method does nothing, so just
  // exit
  return true;
}

