//
// Class Cyclotron
//   Defines the abstract interface for a cyclotron.
//
// Copyright (c) 2007 - 2012, Jianjun Yang and Andreas Adelmann, Paul Scherrer Institut, Villigen PSI, Switzerland
// Copyright (c) 2013 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects"
// and the paper
// "Beam dynamics in high intensity cyclotrons including neighboring bunch effects:
// Model, implementation, and application"
// (https://journals.aps.org/prab/pdf/10.1103/PhysRevSTAB.13.064201)
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
#ifndef CLASSIC_Cyclotron_HH
#define CLASSIC_Cyclotron_HH

#include "AbsBeamline/Component.h"

#include <string>
#include <vector>

class Fieldmap;
class LossDataSink;
class TrimCoil;

struct BfieldData {
    // known from file: field and three theta derivatives
    std::vector<double> bfld_m;   //Bz
    std::vector<double> dbt_m;    //dBz/dtheta
    std::vector<double> dbtt_m;   //d2Bz/dtheta2
    std::vector<double> dbttt_m;  //d3Bz/dtheta3

    // to be calculated in getdiffs: all other derivatives:
    std::vector<double> dbr_m;    // dBz/dr
    std::vector<double> dbrr_m;   // ...
    std::vector<double> dbrrr_m;

    std::vector<double> dbrt_m;
    std::vector<double> dbrrt_m;
    std::vector<double> dbrtt_m;

    // used to get (Br,Btheta,Bz) at any off-plane point
    std::vector<double> f2_m;  // for Bz
    std::vector<double> f3_m;  // for Br
    std::vector<double> g3_m;  // for Btheta

    // Grid-Size
    int nrad_m, ntet_m; // need to be read from inputfile.
    int ntetS_m;        // one more grid line is stored in azimuthal direction
    int ntot_m;         // total grid points number.

    // Mean and Maximas
    double bacc_m, dbtmx_m, dbttmx_m, dbtttmx_m;
};

struct BPositions {
    // these 4 parameters are need to be read from field file.
    double rmin_m, delr_m;
    double tetmin_m, dtet_m;

    // Radii and step width of initial Grid
    std::vector<double> rarr_m;

    // Multiplication factor for magnetic field
    double Bfact_m;
};

class Cyclotron: public Component {

public:
    enum class BFieldType: unsigned short {
        PSIBF,
        CARBONBF,
        ANSYSBF,
        AVFEQBF,
        FFABF,
        BANDRF,
        SYNCHRO
    };

    /// Constructor with given name.
    explicit Cyclotron(const std::string& name);

    Cyclotron();
    Cyclotron(const Cyclotron&);

    virtual ~Cyclotron();

    /// Apply visitor to Cyclotron.
    virtual void accept(BeamlineVisitor&) const;

    /// Get number of slices.
    //  Slices and stepsize used to determine integration step.
    virtual double getSlices() const = 0;

    /// Get stepsize.
    //  Slices and stepsize used to determine integration step.
    virtual double getStepsize() const = 0;

    void setFieldMapFN(const std::string& fmapfn);
    virtual std::string getFieldMapFN() const;

    void setRfFieldMapFN(std::vector<std::string> rffmapfn);
    void setRFFCoeffFN(std::vector<std::string> rff_coeff_fn);
    void setRFVCoeffFN(std::vector<std::string> rfv_coeff_fn);

    void setCyclotronType(const std::string& type);
    const std::string& getCyclotronType() const;

    void setBFieldType();
    BFieldType getBFieldType() const;

    virtual ElementType getType() const;

    virtual void getDimensions(double& zBegin, double& zEnd) const;

    unsigned int getNumberOfTrimcoils() const;

    void setCyclHarm(double h);
    virtual double getCyclHarm() const;

    void setRfPhi(std::vector<double> f);
    virtual std::vector<double> getRfPhi() const;

    void setRfFrequ(std::vector<double> f);
    virtual std::vector<double> getRfFrequ() const;

    void setSymmetry(double symmetry);
    virtual double getSymmetry() const;

    void setRinit(double rinit);
    virtual double getRinit() const;

    void setPRinit(double prinit);
    virtual double getPRinit() const;

    void setPHIinit(double phiinit);
    virtual double getPHIinit() const;

    void setZinit(double zinit);
    virtual double getZinit() const;

    void setPZinit(double zinit);
    virtual double getPZinit() const;

    void setBScale(double bs);
    virtual double getBScale() const;

    void setEScale(std::vector<double> bs);
    virtual std::vector<double> getEScale() const;

    void setTrimCoils(const std::vector<TrimCoil*>& trimcoils);

    void setSuperpose(std::vector<bool> flag);
    virtual std::vector<bool> getSuperpose() const;

    void setMinR(double r);
    virtual double getMinR() const;
    void setMaxR(double r);
    virtual double getMaxR() const;

    void setMinZ(double z);
    virtual double getMinZ() const;
    void setMaxZ(double z);
    virtual double getMaxZ() const;

    void setFMLowE(double e);
    virtual double getFMLowE() const;
    void setFMHighE(double e);
    virtual double getFMHighE() const;

    void setTrimCoilThreshold(double);
    virtual double getTrimCoilThreshold() const;

    void setSpiralFlag(bool spiral_flag);
    virtual bool getSpiralFlag() const;

    virtual bool apply(const size_t& id, const double& t, Vector_t& E, Vector_t& B);

    virtual bool apply(const Vector_t& R, const Vector_t& P, const double& t, Vector_t& E, Vector_t& B);

    virtual void apply(const double& rad, const double& z,
                       const double& tet_rad, double& br,
                       double& bt, double& bz);

    virtual void initialise(PartBunchBase<double, 3>* bunch, double& startField, double& endField);

    virtual void initialise(PartBunchBase<double, 3>* bunch, const double& scaleFactor);

    virtual void finalise();

    virtual bool bends() const;

    virtual double getRmax() const;
    virtual double getRmin() const;

    bool interpolate(const double& rad,
                     const double& tet_rad,
                     double& br,
                     double& bt,
                     double& bz);

    void read(const double& scaleFactor);

    void writeOutputFieldFiles();

private:
    /// Apply trim coils (calculate field contributions) with smooth field transition
    void applyTrimCoil  (const double r, const double z, const double tet_rad, double& br, double& bz);
    /// Apply trim coils (calculate field contributions)
    void applyTrimCoil_m(const double r, const double z, const double tet_rad, double* br, double* bz);


protected:
    void   getdiffs();
    double gutdf5d(double* f, double dx, const int kor, const int krl, const int lpr);

    void   initR(double rmin, double dr, int nrad);

    void   getFieldFromFile_Ring(const double& scaleFactor);
    void   getFieldFromFile_Carbon(const double& scaleFactor);
    void   getFieldFromFile_CYCIAE(const double& scaleFactor);
    void   getFieldFromFile_AVFEQ(const double& scaleFactor);
    void   getFieldFromFile_FFA(const double& scaleFactor);
    void   getFieldFromFile_BandRF(const double& scaleFactor);
    void   getFieldFromFile_Synchrocyclotron(const double& scaleFactor);

    inline int idx(int irad, int ktet) {return (ktet + Bfield_m.ntetS_m * irad);}


private:
    BFieldType fieldType_m;

    std::string fmapfn_m; /**< Stores the filename of the B-fieldmap*/
    std::vector<double> rffrequ_m;
    std::vector< std::vector<double> > rffc_m;
    std::vector<double> rfvrequ_m;
    std::vector< std::vector<double> > rfvc_m;
    std::vector<double> rfphi_m;
    std::vector<double> escale_m;  /**< A scale factor for the E-field*/
    std::vector<bool> superpose_m; /**< A flag for superpose electric fields*/

    double symmetry_m;

    double rinit_m;
    double prinit_m;
    double phiinit_m;
    double zinit_m;
    double pzinit_m;

    bool spiralFlag_m;
    double trimCoilThreshold_m; /**< B-field threshold for applying trim coil*/

    std::string typeName_m; /**< Name of the TYPE parameter in cyclotron*/

    double harm_m;

    double bscale_m; /**< A scale factor for the B-field*/

    std::vector<TrimCoil*> trimcoils_m; /**< Trim coils*/

    double minr_m;
    double maxr_m;
    double minz_m;
    double maxz_m;

    double fmLowE_m;
    double fmHighE_m;

    // Not implemented.
    void operator=(const Cyclotron &) = delete;

    // RF field map handler
    //    Fieldmap *RFfield;
    std::vector<Fieldmap*> RFfields_m;
    std::vector<std::string> RFfilename_m;
    std::vector<std::string> RFFCoeff_fn_m;
    std::vector<std::string> RFVCoeff_fn_m;

    std::unique_ptr<LossDataSink> lossDs_m; /**< Handling for store the particle out of region*/

    // Necessary for quick and dirty phase output -DW
    int waitingGap_m = 1;

protected:
    // object of Matrices including magnetic field map and its derivates
    BfieldData Bfield_m;

    // object of parameters about the map grid
    BPositions BP_m;
};

#endif // CLASSIC_Cyclotron_HH
