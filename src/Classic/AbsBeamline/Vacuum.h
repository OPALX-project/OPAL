//
// Class Vacuum
//   Defines the abstract interface for vacuum.
//
// Copyright (c) 2018 - 2021, Pedro Calvo, CIEMAT, Spain
// All rights reserved.
//
// Implemented as part of the PhD thesis
// "Optimizing the radioisotope production of the novel AMIT
// superconducting weak focusing cyclotron"
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
#ifndef CLASSIC_Vacuum_HH
#define CLASSIC_Vacuum_HH

#include "AbsBeamline/Component.h"
#include "AbsBeamline/Cyclotron.h"

#include <boost/bimap.hpp>

#include <memory>
#include <string>
#include <vector>

class BeamlineVisitor;
class Cyclotron;

struct PFieldData {
    std::vector<double> pfld_m; // Pressure field from file
    int nrad_m, ntet_m;         // Grid-Size, from inputfile.
    int ntetS_m;                // extra grid line is stored in azimuthal direction
    int ntot_m;                 // total grid points number
};

struct PPositions {
    // these 4 parameters are need to be read from field file.
    double rmin_m, delr_m;
    double tetmin_m, dtet_m;

    // Radii and step width of initial Grid
    std::vector<double> rarr_m;

    double Pfact_m;
};

enum class ResidualGas: short {
    NOGAS = -1,
    AIR   = 0,
    H2    = 1
};

class Vacuum: public Component {

public:
    /// Constructor with given name.
    explicit Vacuum(const std::string& name);

    Vacuum();
    Vacuum(const Vacuum& rhs);
    virtual ~Vacuum();

    /// Apply visitor to Vacuum.
    virtual void accept(BeamlineVisitor&) const override;

    virtual bool apply(const size_t& i, const double& t, Vector_t& E, Vector_t& B) override;

    virtual bool applyToReferenceParticle(const Vector_t& R, const Vector_t& P,
                                          const double& t, Vector_t& E, Vector_t& B) override;

    virtual bool checkVacuum(PartBunchBase<double, 3>* bunch, Cyclotron* cycl);

    virtual void initialise(PartBunchBase<double, 3>* bunch,
                            double& startField, double& endField) override;

    virtual void initialise(PartBunchBase<double, 3>* bunch);

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual void goOnline(const double& kineticEnergy) override;

    virtual void goOffline() override;

    virtual ElementType getType() const override;

    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    bool checkPoint(const Vector_t& R);

    double checkPressure(const Vector_t& R);

    void setResidualGas(std::string gas);
    ResidualGas getResidualGas() const;
    std::string getResidualGasName();

    void setPressure(double pressure);
    double getPressure() const;

    void setPressureMapFN(std::string pmapfn);
    std::string getPressureMapFN() const;

    void setPScale(double ps);
    double getPScale() const;

    void setTemperature(double temperature);
    double getTemperature() const;

    void setStop(bool stopflag);
    bool getStop() const;


protected:
    void initR(double rmin, double dr, int nrad);

    void getPressureFromFile();

    inline int idx(int irad, int ktet) {return (ktet + PField_m.ntetS_m * irad);}


private:
    // Not implemented.
    void operator=(const Vacuum&);

    void updateParticleAttributes();

    void print();

    ///@{ parameters for Vacuum
    ResidualGas gas_m;    /// Type of gas for residual vacuum
    double pressure_m;    /// mbar
    std::string pmapfn_m; /// stores the filename of the pressure map
    double pscale_m;      /// a scale factor for the P-field
    double temperature_m; /// K
    bool stop_m;          /// Flag if particles should be stripped or stopped
    ///@}

    ///@{ size limits took from cyclotron
    double minr_m;   /// mm
    double maxr_m;   /// mm
    double minz_m;   /// mm
    double maxz_m;   /// mm
    ///@}

    ParticleMatterInteractionHandler* parmatint_m;

    static const boost::bimap<ResidualGas, std::string> bmResidualGasString_s;


protected:
    // object of Matrices including pressure field map and its derivates
    PFieldData PField_m;

    // object of parameters about the map grid
    PPositions PP_m;
};


inline
bool Vacuum::bends() const {
    return false;
}

#endif // CLASSIC_Vacuum_HH
