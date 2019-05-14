#ifndef AMR_OBECT_H
#define AMR_OBECT_H

#include "Index/NDIndex.h"

#include "Algorithms/PBunchDefs.h"
#include "Utilities/Util.h"

// #include "Algorithms/AmrPartBunch.h"

// class AmrPartBunch;

/**
 * The AMR interface to OPAL. A new AMR library needs
 * to inherit from this class in order to work properly
 * with OPAL. Among other things it specifies the refinement
 * strategies.
 */
class AmrObject {
    
public:
    // FIXME Why not using typedef of PartBunchBase
    typedef std::pair<Vector_t, Vector_t> VectorPair_t;
//     using VectorPair_t = typename AmrPartBunch::VectorPair_t;
    
    /// Methods for tagging cells for refinement
    enum TaggingCriteria {
        CHARGE_DENSITY = 0, // default
        POTENTIAL,
        EFIELD,
        MOMENTA,
        MIN_NUM_PARTICLES,       ///< min. #particles per cell
        MAX_NUM_PARTICLES        ///< max. #particles per cell
    };
    
    
    /*!
     * This data structure is only used for creating an object
     * via the static member function AmrBoxLib::create()
     * that is called in FieldSolver::initAmrObject_m
     */
    struct AmrInfo {
        int grid[3];        ///< Number of grid points in x-, y- and z-direction
        int maxgrid[3];     ///< Maximum grid size in x-, y- and z-direction
        int bf[3];          ///< Grid blocking factor in x-, y- and z-direction
        int maxlevel;       ///< Maximum level for AMR (0: single-level)
        int refratio[3];    ///< Mesh refinement ratio in x-, y- and z-direction
    };
    
public:
    
    AmrObject();
    
    AmrObject(TaggingCriteria tagging,
              double scaling,
              double chargedensity);
    
    virtual ~AmrObject();
    
    /*!
     * Collect information about grid load balancing.
     * @param gridPtsPerCore is filled.
     * @param gridsPerLevel is filled
     */
    virtual void getGridStatistics(std::map<int, long>& gridPtsPerCore,
                                   std::vector<int>& gridsPerLevel) const = 0;
    
    /*!
     * Setup all fine levels after object creation.
     */
    virtual void initFineLevels() = 0; 
    
    /*!
     * Update of mesh according to chosen refinement strategy.
     * @param time of regrid
     */
    virtual void regrid(double time) = 0;
    
    /*!
     * Choose a new tagging strategy.
     * Is used in src/Structure/FieldSolver.cpp
     * @param tagging strategy
     */
    void setTagging(TaggingCriteria tagging);
    
    /*!
     * Choose a new tagging strategy (string version).
     * Is used in src/Structure/FieldSolver.cpp
     * @param tagging strategy
     */
    void setTagging(std::string tagging);
    
    /*!
     * Scaling factor for tagging.
     * It is used with POTENTIAL and EFIELD
     * @param scaling factor in [0, 1]
     */
    void setScalingFactor(double scaling);
    
    /*!
     * Charge density for tagging with CHARGE_DENSITY
     * @param chargedensity >= 0.0 (e.g. 1e-14)
     */
    void setChargeDensity(double chargedensity);
    
    /*!
     * Maximum number of particles per cell for tagging
     * @param maxNumPart is upper bound for a cell to be marked
     * for refinement
     */
    void setMaxNumParticles(size_t maxNumPart);
    
    /*!
     * Minimum number of particles per cell for tagging
     * @param minNumPart is lower bound for a cell to be marked
     * for refinement
     */
    void setMinNumParticles(size_t minNumPart);
    
    /* Methods that are needed by the
     * bunch
     */
    virtual VectorPair_t getEExtrema() = 0;
    
    virtual double getRho(int x, int y, int z) = 0;
    
    virtual void computeSelfFields() = 0;
    
    virtual void computeSelfFields(int b) = 0;
    
    virtual void computeSelfFields_cycl(double gamma) = 0;
    
    virtual void computeSelfFields_cycl(int b) = 0;
    
    virtual void updateMesh() = 0;
    
    virtual Vektor<int, 3> getBaseLevelGridPoints() const = 0;
    
    virtual const int& maxLevel() const = 0;
    virtual const int& finestLevel() const = 0;
    
    /*!
     * @returns the time of the simulation
     */
    virtual double getT() const = 0;
    
    /*!
     * Rebalance the grids among the
     * cores
     */
    virtual void redistributeGrids(int how) { }
    
    /*!
     * Used in AmrPartBunch to check if we need to refine
     * first.
     * @returns true fine grids are initialized
     */
    const bool& isRefined() const;
    
    /*!
     * Used in Fieldsolver in order to convert a number that
     * specifies the tagging to the corresponding string. Check
     * enum TaggingCriteria for ordering.
     * @param number of tagging
     */
    static std::string enum2string(int number);

protected:

    TaggingCriteria tagging_m;  ///< Tagging strategy
    
    double scaling_m;           ///< Scaling factor for tagging [0, 1]
                                // (POTENTIAL, EFIELD)
    double chargedensity_m;     ///< Tagging value for CHARGE_DENSITY
    
    size_t maxNumPart_m;        ///< Tagging value for MAX_NUM_PARTICLES
    
    size_t minNumPart_m;        ///< Tagging value for MIN_NUM_PARTICLES
    
    bool refined_m;             ///< Only set to true in AmrObject::initFineLevels()
    
    /// timer for selfField calculation (used in concrete AmrObject classes)
    IpplTimings::TimerRef amrSolveTimer_m;
    IpplTimings::TimerRef amrRegridTimer_m;
};

#endif
