#ifndef NO_PARTICLE_CACHING_POLICY
#define NO_PARTICLE_CACHING_POLICY

/*
 * 
 * Empty caching strategy that doesn't cache anything
 * 
 */

template <class T, unsigned Dim, class Mesh, class CachingPolicy> class ParticleSpatialLayout;

//basic policy that doesn't cache any particles
template<class T, unsigned Dim, class Mesh>
class NoParticleCachingPolicy {
public:
template<class C>
	void updateCacheInformation(
		ParticleSpatialLayout<T, Dim, Mesh, C > &/*PLayout*/
		)
	{
		//don't do anything...
	}
template<class C>
	void updateGhostParticles(
		IpplParticleBase< ParticleSpatialLayout<T,Dim,Mesh,C > > &/*PData*/,
		ParticleSpatialLayout<T, Dim, Mesh, C > &/*PLayout*/
		)
	{
		//don't do anything...
	}
protected:
	~NoParticleCachingPolicy() {}
};

#endif
