#ifndef OPAL_SECONDARY_EMISSION_PHYSICS_HH
#define OPAL_SECONDARY_EMISSION_PHYSICS_HH

#include <hdf5.h>
#include <limits>

#include "Algorithms/PBunchDefs.h"
#include "Algorithms/PartBunchBase.h"

extern Inform *gmsg;

namespace myeps {
    const double EPS = std::numeric_limits<double>::epsilon();
    const double FPMIN = std::numeric_limits<double>::min() / EPS;
}

class SecondaryEmissionPhysics {

public:

    SecondaryEmissionPhysics();

    ~SecondaryEmissionPhysics();

    void nSec(const double &incEnergy,
              const double &cosTheta,
              const int &matNumber,
              int &seNum,
              int &seType,
              const double &incQ,
              const Vector_t &TriNorm,
              const Vector_t &inteCoords,
              const Vector_t &localX,
              PartBunchBase<double, 3> *itsBunch,
              double &seyNum, const
              double &ppVw,
              const double &vVThermal,
              const bool nEmissionMode);

    void nSec(const double &incEnergy,
              const double &cosTheta,
              int &seNum,
              int &seType,
              const double &incQ,
              const Vector_t &TriNorm,
              const Vector_t &inteCoords,
              const Vector_t &localX,
              PartBunchBase<double, 3> *itsBunch,
              double &seyNum,
              const double &ppVw,
              const double &vSeyZero,
              const double &vEzero,
              const double &vSeyMax,
              const double &vEmax,
              const double &vKenergy,
              const double &vKtheta,
              const double &vVThermal,
              const bool nEmissionMode);

    double deltae_m;
    double deltar_m;
    double deltats_m;

    private:

    /**
     * @param TPnSec_m is the timmer of secondary emission module.
     */

    IpplTimings::TimerRef TPnSec_m;

    /**
     * @param seAlpha_m is the emitted angular spectrum parameter, i.e., the 1st parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seAlpha_m;

    /*==============================================================================================
      parameters for backscattered electrons
      ==============================================================================================*/

    /**
     * @param sePScat_m is the 2nd parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double sePScat_m;
    /**
     * @param sePScatPeak_m is the 3rd parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double sePScatPeak_m;
    /**
     * @param seEScatPeak_m is the 4th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seEScatPeak_m;
    /**
     * @param seW_m is the 5th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seW_m;
    /**
     * @param seP_m is the 6th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seP_m;
    /**
     * @param seDelta_m is the 7th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seDelta_m;
    /**
     * @param seEOne_m is the 8th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seEOne_m;
    /**
     * @param seETwo_m is the 9th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seETwo_m;


    /*==============================================================================================
      parameters for rediffused electrons
      ==============================================================================================*/

    /**
     * @param sePRed_m is the 10th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double sePRed_m;
    /**
     * @param seERed_m is the 11th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seERed_m;
    /**
     * @param seR_m is the 12th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seR_m;
    /**
     * @param seQ_m is the 13th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seQ_m;
    /**
     * @param seROne_m is the 14th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seROne_m;
    /**
     * @param seRTwo_m is the 15th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seRTwo_m;


    /*==============================================================================================
      parameters for true secondary electrons
      ==============================================================================================*/

    /**
     * @param seYPeakTS_m is the 16th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seYPeakTS_m;
    /**
     * @param seEPeakTS_m is the 17th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seEPeakTS_m;
    /**
     * @param seSTS_m is the 18th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seSTS_m;
    /**
     * @param seTOneTS_m is the 19th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seTOneTS_m;
    /**
     * @param seTTwoTS_m is the 20th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seTTwoTS_m;
    /**
     * @param seTThreeTS_m is the 21st parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seTThreeTS_m;
    /**
     * @param seTFourTS_m is the 22nd parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seTFourTS_m;
    /**
     * @param seEPeakTot_m is the 23rd parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seEPeakTot_m;
    /**
     * @param seYPeakTot_m is the 24th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seYPeakTot_m;

    /*==============================================================================================
      parameters in he TABLE II of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
      ==============================================================================================*/
    double sePn_m[10];
    double seEpsn_m[10];


    void calcEmiNum(const double incEnergy,
                    const double cosTheta,
                    int &seNum,
                    const double &vSeyZero,
                    const double &vEzero,
                    const double &vSeyMax,
                    const double &vEmax,
                    const double &vKenergy,
                    const double &vKtheta,
                    double &seyNum);

    void calcEmiNum(const double incEnergy, const double cosTheta, const double *prob, int &seNum);

    double calcDeltats(const double incEnergy, const double cosTheta);

    double calcDeltar(const double incEnergy, const double cosTheta);

    double calcDeltae(const double incEnergy, const double cosTheta);

    double calcProb(const double incEnergy, const double cosTheta, double *prob);

    void setSeMaterial(int material_num);

    /*==========================================
     *http://www.taygeta.com/random/gaussian.html
     *return a gaussian distributed random number
     *==========================================*/

    double gaussRand();

    double gammp(const double a, const double x);

    double gser(const double a, const double x);
    double gcf(const double a, const double x);
    double gammpapprox(double a, double x, int psig);
    double gammln(const double xx);

    double invgammp(double p, double a);
    double min_Gamma(double x, double y);
    double max_Gamma(double x, double y);

    double betai(const double x, const double a, const double b);
    double betacf(const double a, const double b, const double x);
    double betaiapprox(double a, double b, double x);
    double invbetai(double p, double a, double b);
    void coordConverter (const Vector_t &/*TriNormal*/, Vector_t &/*x*/) {  }

};

inline
double SecondaryEmissionPhysics::min_Gamma(double x, double y) {
    if (x > y)
        return y;
    else
        return x;
}

inline
double SecondaryEmissionPhysics::max_Gamma(double x, double y) {
    if (x < y)
        return y;
    else
        return x;
}

/*  ==========================================================================*/

#endif //OPAL_SECONDARY_EMISSION_PHYSICS_HH