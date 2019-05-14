#ifndef AMR_PART_BUNCH_H
#define AMR_PART_BUNCH_H

#include "Algorithms/PartBunchBase.h"
#include "Amr/AmrObject.h"

class AmrPartBunch : public PartBunchBase<double, 3>
{
public:
    typedef AmrParticle_t pbase_t;
    
public:

    AmrPartBunch(const PartData *ref);

    AmrPartBunch(const PartData *ref, pbase_t* pbase_p);

    ~AmrPartBunch();

    pbase_t *getAmrParticleBase();

    const pbase_t *getAmrParticleBase() const;

    void initialize(FieldLayout_t *fLayout);
    
    // does actually another repartition
    void do_binaryRepart();
    
    Vector_t get_hr() const;
    
    void set_meshEnlargement(double dh);
    
    VectorPair_t getEExtrema();
    
    double getRho(int x, int y, int z);
    
    FieldLayout_t &getFieldLayout();
    
    
    void boundp();
    
    void computeSelfFields();
    
    void computeSelfFields(int bin);
    
    void computeSelfFields_cycl(double gamma);
    
    void computeSelfFields_cycl(int bin);
    
    void setSolver(FieldSolver *fs) {
        PartBunchBase<double, 3>::setSolver(fs);
        this->amrobj_mp = fs->getAmrObject();
    }
    
    virtual void setBinCharge(int bin, double q) { };
    virtual void setBinCharge(int bin) { };
    
    /*
     * AmrPartBunch only
     */
    
    const AmrObject* getAmrObject() const {
        return this->amrobj_mp;
    }
    
    PoissonSolver *getFieldSolver() {
        return fs_m->solver_m;
    }
    
    const PoissonSolver *getFieldSolver() const {
        return fs_m->solver_m;
    }
    
    void setBaseLevelMeshSpacing(const Vector_t& hr) {
        for (int i = 0; i < 3; ++i)
            hr_m[i] = hr[i];
    }
    
    void gatherLevelStatistics();
    
    /*!
     * Only a valid call of root core (core 0)
     * @param l is the level
     */
    const size_t& getLevelStatistics(int l) const;
    
    
    /*!
     * Update the Lorentz factor before every domainMapping in order
     * to have the correct boosted frame for the particle redistribution,
     * regrid and computation of the self-field forces
     * @param bin is the energy bin
     */
    void updateLorentzFactor(int bin=0);
    
    /*!
     * Update the Lorentz factor before every domainMapping in order
     * to have the correct boosted frame for the particle redistribution,
     * regrid and computation of the self-field forces. This function is
     * only used in the single bunch case.
     * @param gamma is the Lorentz factor
     */
    void updateLorentzFactor(double gamma);
    
    //FIXME BCs
    void setBCAllPeriodic() {}
    void setBCAllOpen() {}
    void setBCForDCBeam() {}
    
    
private:
    void updateFieldContainers_m();
    
    void updateDomainLength(Vektor<int, 3>& grid);
    
    void updateFields(const Vector_t& hr, const Vector_t& origin);
    
private:
    
    /* pointer to AMR object that is part
     * of solver_m (AmrPoissonSolver) in src/Structure/FieldSolver.h
     */
    AmrObject *amrobj_mp;
    pbase_t *amrpbase_mp;
    
    /* We need this due to H5PartWrapper etc, but it's always nullptr.
     * Thus, don't use it.
     */
    FieldLayout_t* fieldlayout_m;
    
    std::unique_ptr<size_t[]> globalPartPerLevel_m;
};

#endif
