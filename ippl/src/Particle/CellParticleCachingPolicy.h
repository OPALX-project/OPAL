#ifndef CELL_PARTICLE_CACHING_POLICY
#define CELL_PARTICLE_CACHING_POLICY

/*
 *
 * The Cell caching layout ensures that each node has all ghost particles
 * for each external particle that is inside a neighboring cell.
 *
 */

#include <Particle/BoxParticleCachingPolicy.h>

template<class T, unsigned Dim, class Mesh>
class CellParticleCachingPolicy : private BoxParticleCachingPolicy<T,Dim,Mesh> {
public:
	CellParticleCachingPolicy()
	{
		std::fill(cells, cells+Dim, 0);
	}

	void setCacheCellRange(int d, int length)
	{
		cells[d] = length;
	}

	void setAllCacheCellRanges(int length)
	{
		std::fill(cells, cells+Dim, length);
	}

	template<class C>
	void updateCacheInformation(
		ParticleSpatialLayout<T, Dim, Mesh, C > &PLayout
		)
	{
		for(unsigned int d = 0;d<Dim;++d)
			BoxParticleCachingPolicy<T,Dim,Mesh>::setCacheDimension(d, cells[d]*PLayout.getLayout().getMesh().get_meshSpacing(d));

		BoxParticleCachingPolicy<T,Dim,Mesh>:: updateCacheInformation(PLayout);
	}

	template<class C>
	void updateGhostParticles(
		IpplParticleBase< ParticleSpatialLayout<T,Dim,Mesh,C > > &PData,
		ParticleSpatialLayout<T, Dim, Mesh, C > &PLayout
		)
	{
		for(unsigned int d = 0;d<Dim;++d)
			BoxParticleCachingPolicy<T,Dim,Mesh>::setCacheDimension(d, cells[d]*PLayout.getLayout().getMesh().get_meshSpacing(d));

		BoxParticleCachingPolicy<T,Dim,Mesh>::updateGhostParticles(PData, PLayout);
	}
protected:
	~CellParticleCachingPolicy() {}
private:
	int cells[Dim];
};

#endif