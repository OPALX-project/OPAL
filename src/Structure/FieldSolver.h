//
// Class FieldSolver
//   The class for the OPAL FIELDSOLVER command.
//   A FieldSolver definition is used by most physics commands to define the
//   particle charge and the reference momentum, together with some other data.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
//
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
#ifndef OPAL_FieldSolver_HH
#define OPAL_FieldSolver_HH

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Solvers/PoissonSolver.h"
#ifdef ENABLE_AMR
    #include "Amr/AmrObject.h"
    #include "Solvers/AmrPoissonSolver.h"
    #include <memory>
#endif

#include <string>

template <class T, unsigned Dim>
class PartBunchBase;

enum class FieldSolverType: short {
    NONE = -1,
    FFT,
    FFTBOX,
    SAAMG,
    P3M,
    FMG,
    ML,
    AMRMG,
    HYPRE,
    HPGMG
};


class FieldSolver: public Definition {

public:
    /// Exemplar constructor.
    FieldSolver();

    virtual ~FieldSolver();

    /// Make clone.
    virtual FieldSolver* clone(const std::string& name);

    /// Find named FieldSolver.
    static FieldSolver* find(const std::string& name);

    std::string getType();

    /// Return meshsize
    double getMX() const;

    /// Return meshsize
    double getMY() const;

    /// Return meshsize
    double getMT() const;

    /// Store emittance for mode 1.
    void setMX(double);

    /// Store emittance for mode 2.
    void setMY(double);

    /// Store emittance for mode 3.
    void setMT(double);

    /// Update the field solver data.
    virtual void update();

    /// Execute (init) the field solver data.
    virtual void execute();

    void initCartesianFields();

    void initSolver(PartBunchBase<double, 3>* b);

    bool hasValidSolver();

    void setFieldSolverType();
    FieldSolverType getFieldSolverType() const;

    inline Layout_t &getParticleLayout() { return* PL_m; }

    FieldLayout_t *getFieldLayout() { return FL_m; }

    Inform& printInfo(Inform& os) const;

    unsigned int getInteractionRadius() {return (unsigned int) rpp_m; }

    bool hasPeriodicZ();

    bool isAmrSolverType() const;

#ifdef ENABLE_AMR
    AmrObject *getAmrObject() {
        return itsAmrObject_mp.get();
    }
    
    const AmrObject *getAmrObject() const {
        return itsAmrObject_mp.get();
    }
#endif

    /// the actual solver, should be a base object
    PoissonSolver* solver_m;

private:
#ifdef ENABLE_AMR

    std::string getTagging_m() const;

    void initAmrObject_m();
    
    void initAmrSolver_m();
    
    std::unique_ptr<AmrObject> itsAmrObject_mp;
#endif

    // Not implemented.
    FieldSolver(const FieldSolver&);
    void operator=(const FieldSolver&);

    // Clone constructor.
    FieldSolver(const std::string& name, FieldSolver* parent);

    /// The cartesian mesh
    Mesh_t* mesh_m;

    /// The field layout f
    FieldLayout_t* FL_m;

    /// The particle layout
    std::unique_ptr<Layout_t> PL_m;

    /// all the particles are here ...
    PartBunchBase<double, 3>* itsBunch_m;

    std::string fsName_m;
    FieldSolverType fsType_m;

    double rpp_m;
};

inline
FieldSolverType FieldSolver::getFieldSolverType() const {
    return fsType_m;
}

inline
Inform& operator<<(Inform& os, const FieldSolver& fs) {
    return fs.printInfo(os);
}

#endif // OPAL_FieldSolver_HH
