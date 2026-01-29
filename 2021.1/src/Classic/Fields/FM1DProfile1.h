#ifndef CLASSIC_FIELDMAP1DPROFILE1_HH
#define CLASSIC_FIELDMAP1DPROFILE1_HH

#include "Fields/Fieldmap.h"

/*
 * Class FM1DProfile.
 * --------------------------------------------------------------------------
 * Field definition for 1D representation of bending magnet.
 *
 * Class FM1DProfile1 defines a 1D field map for us in bending magnets.
 */

class FM1DProfile1: public Fieldmap {

public:

    virtual bool getFieldDerivative(const Vector_t &X,
                                    Vector_t &E,
                                    Vector_t &B,
                                    const DiffDirection &dir) const;
    virtual void get1DProfile1EntranceParam(double &entranceParameter1,
                                           double &entranceParameter2,
                                           double &entranceParameter3);
    virtual void get1DProfile1ExitParam(double &exitParameter1,
                                       double &exitParameter2,
                                       double &exitParameter3);
    virtual double getFieldGap();
    virtual void getFieldDimensions(double &zBegin,
                                    double &zEnd) const;
    virtual void getFieldDimensions(double &xIni,
                                    double &xFinal,
                                    double &yIni,
                                    double &yFinal,
                                    double &zIni,
                                    double &zFinal) const;
    virtual bool getFieldstrength(const Vector_t &X,
                                  Vector_t &strength,
                                  Vector_t &info) const;
    virtual double getFrequency() const;
    virtual void getInfo(Inform *);
    virtual void setFrequency(double freq);
    virtual void get1DProfile1EngeCoeffs(std::vector<double> &engeCoeffsEntry,
                                         std::vector<double> &engeCoeffsExit);
    virtual void swap();

    virtual void setFieldGap(double gap);

private:

    /// Constructor with field map file name.
    FM1DProfile1(std::string Filename);

    virtual ~FM1DProfile1();

    virtual void freeMap();
    virtual void readMap();

    double computeEntranceFringe(double z) const;
    double computeExitFringe(double z) const;
    double computeFringe(const std::vector<double> &coefs, double z) const;
    /*
     * Entrance and exit position parameters. These are read in from the input
     * input file. Ultimately they are used to determine the origin of the
     * field Enge function and the extent of the field map. However, how they
     * are used to do this depends on how the bend using the map is setup in the
     * OPAL input file. So, we use generic terms to start.
     */
    double entranceParameter1_m;
    double entranceParameter2_m;
    double entranceParameter3_m;
    double exitParameter1_m;
    double exitParameter2_m;
    double exitParameter3_m;

    /// Enge coefficients for map entry and exit regions.
     std::vector<double> engeCoeffsEntry_m;
     std::vector<double> engeCoeffsExit_m;

     int polyOrderEntry_m;           /// Enge function order for entry region.
     int polyOrderExit_m;            /// Enge function order for entry region.
     double gapHeight_m;             /// Full gap height of field map.

     double sBegin_m;                /// Start of field map in s coordinates (m).
     double sEnd_m;                  /// End of field map in s coordinates (m).

    friend class Fieldmap;
};

#endif