#include "Structure/SecondaryEmissionPhysics.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>
#include "Utilities/Options.h"
#include "Utility/IpplTimings.h"

using namespace myeps;
using namespace Physics;
SecondaryEmissionPhysics::SecondaryEmissionPhysics() {

    TPnSec_m = IpplTimings::getTimer("Secondary emission");
}
/**
 * Destructor.
 *
 * Delete the previous defined member arrays
 */
SecondaryEmissionPhysics::~SecondaryEmissionPhysics() {

}

void SecondaryEmissionPhysics::nSec(const double &incEnergy,
                                    const double &cosTheta,
                                    const int &matNumber,
                                    int &seNum,
                                    int &seType,
                                    const double &incQ,
                                    const Vector_t &TriNorm,
                                    const Vector_t &inteCoords,
                                    const Vector_t &localX,
                                    PartBunchBase<double, 3> *itsBunch,
                                    double &seyNum,
                                    const double &ppVw,
                                    const double &vVThermal,
                                    const bool nEmissionMode) {

    IpplTimings::startTimer(TPnSec_m);
    double prob[11] = {0};
    std::vector<Vector_t> se_P;
    setSeMaterial(matNumber);//set material based parameters
    seyNum=calcProb(incEnergy, cosTheta, prob);//calculate probability
    calcEmiNum(incEnergy, cosTheta, prob, seNum);//calculate emitted number
    PAssert_LT(seNum, 11);//debug
    PAssert_GE(seNum, 0);//debug
    double Eemit[10];
    double emiTheta[10];
    double emiPhi[10];
    Vector_t interCoords_l = inteCoords;
    Vector_t TriNorm_l = TriNorm;
    double incQ_l = incQ;


    /*===========================Definitions for benchmark===================================*/
    double vw=ppVw; //1.6*1e-19*1200/9.10938188*1e-31/(2*3.1415926*2.0*1e8)/0.03;//benchmark
    double vt=vVThermal;//7.268929821*1e5;//1.5eV//benchmark
    double f_max=vw/vt*exp(-0.5);//benchmark
    double test_a=vt/vw;//benchmark
    double test_asq=test_a*test_a;//benchmark
    /*---------------------------------------------------------------------------------------*/
    if( Options::ppdebug ) {

    } else {

        if(seNum != 0) {

            for(int i = 0; i < seNum; i++) {

                double tmp1 = IpplRandom();
                double tmp2 = IpplRandom();
                double temp = 1.0 / (1.0 + seAlpha_m);
                emiTheta[i] = acos(pow(tmp1, temp));
                emiPhi[i] = Physics::two_pi * tmp2;

            }
        }
    }


    if(seNum == 0) {

        // The incident particle will be marked for deletion
        if (!nEmissionMode) {
            if( Options::ppdebug ) {
                /*=========(Velocity with Maxwellian Distribution For Parallel Plate Benchmark)========*/
                double test_s=1;
                double f_x=0;
                double test_x=0;
                while (test_s>f_x) {
                    test_s=IpplRandom();
                    test_s*=f_max;
                    test_x=IpplRandom();
                    test_x*=10*test_a;// range for normalized emission speed(0,10*test_a);
                    f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
                }
                double v_emi=test_x*vw;
                Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;
                // cout<<"Single Eemit[0]: "<<Eemit[0]<<endl;
                /*---------------------End Of Maxwellian Distribution(For Benchmark)--------------------*/
            } else {

                // For absorption case in constant simulation particle mode, just use true secondary emission with 1 particle energy distribution for Furman Pivi model.
                double u2 = IpplRandom();

                double temp = incEnergy / seEpsn_m[0];
                double p0 = gammp(sePn_m[0], temp);
                temp = p0 * u2;
                Eemit[0] = seEpsn_m[0] * invgammp(temp, sePn_m[0]) ;


            }

        }


    } else if(seNum == 1) {

        if( Options::ppdebug ) {
            /*=========(Velocity with Maxwellian Distribution For Parallel Plate Benchmark)========*/
            double test_s=1;
            double f_x=0;
            double test_x=0;
            while (test_s>f_x) {
                test_s=IpplRandom();
                test_s*=f_max;
                test_x=IpplRandom();
                test_x*=10*test_a;//range for normalized emission speed(0,10*test_a);
                f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
            }
            double v_emi=test_x*vw;
            Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;
            //cout<<"Single Eemit[0]: "<<Eemit[0]<<endl;
            /*---------------------End Of Maxwellian Distribution(For Benchmark)--------------------*/
        } else {

            //double delta_e = calcDeltae(incEnergy, cosTheta);
            //double delta_r = calcDeltar(incEnergy, cosTheta);
            double tmp = prob[1] + deltae_m + deltar_m;
            //double ae = delta_e / tmp;
            double ae = deltae_m / tmp;
            //double ar = delta_r / tmp;
            double ar = deltar_m / tmp;
            double a_re = ae + ar;
            double urand = IpplRandom();


            if(urand < ae) {
                int t_count = 0;
                do {
                    Eemit[0] = incEnergy - seDelta_m * fabs(gaussRand()) ;
                    t_count++;
                } while(Eemit[0] < 0&&t_count<200);
                if(Eemit[0]<0)// if the above do - while loops over 200 times, the loop will break out, and Eemit will be its mean value, i.e., incident energy.
                    Eemit[0]=incEnergy;
                seType = 0;

            } else if(urand >= ae && urand < a_re) {
                double u1 = IpplRandom();

                double powArg = 1.0 / (1.0 + seQ_m);
                Eemit[0] = incEnergy * pow(u1, powArg);
                seType = 1;


            } else {

                double u2 = IpplRandom();

                double temp = incEnergy / seEpsn_m[0];
                double p0 = gammp(sePn_m[0], temp);
                temp = p0 * u2;
                Eemit[0] = seEpsn_m[0] * invgammp(temp, sePn_m[0]) ;
                seType = 2;

            }
        }

    } else {
        seType = 2;

        if( Options::ppdebug ) {
            /*==========(Velocity with Maxwellian Distribution For Parallel Plate Benchmark)========*/
            if (!nEmissionMode) {
                /*double Eemit_mean = 0.0;
                  for(int i = 0; i < seNum; i++) {
                  Eemit_mean += Eemit[i];
                  }*/
                //Eemit[0] = Eemit_mean/seNum;
                double test_s=1.0;
                double f_x=0;
                double test_x=0;
                while (test_s>f_x) {
                    test_s=IpplRandom();
                    test_s*=f_max;
                    test_x=IpplRandom();
                    test_x*=10*test_a;//range for normalized emission speed(0,10*test_a);
                    f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
                }
                double v_emi=test_x*vw;
                Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;


            } else {
                for(int i = 0; i < seNum; i++) {
                    double test_s=1.0;
                    double f_x=0;
                    double test_x=0;
                    while (test_s>f_x) {
                        test_s=IpplRandom();
                        test_s*=f_max;
                        test_x=IpplRandom();
                        test_x*=10*test_a;//range for normalized emission speed(0,10*test_a);
                        f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
                    }
                    double v_emi=test_x*vw;
                    Eemit[i]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;

                }


            }
            /*---------------------End Of Maxwellian Distribution(For Benchmark)--------------------*/
        } else {
            /*=======================3D velocity for Furman-Pivi's Model====================*/
            double x0 = incEnergy / seEpsn_m[seNum-1];
            double parg = seNum * sePn_m[seNum-1];
            double p0 =  gammp(parg, x0);
            double sin2theta_n[seNum];
            double cos2theta_n[seNum];
            double rand_y =  IpplRandom();
            double invarg = rand_y * p0;
            double y2 = invgammp(invarg, parg);
            double y2_n[seNum];
            double multisin = 1.0;

            if (!nEmissionMode) {// only emit 1 particle

                double Eemisum = 0.0;
                for(int i = 0; i < seNum - 1; i++) {
                    double mu = sePn_m[seNum-1] * (seNum  - i);
                    double nu =  sePn_m[seNum-1];
                    double rand_n = IpplRandom();
                    sin2theta_n[i] = invbetai(rand_n, mu, nu);
                    cos2theta_n[i] = 1.0 - sin2theta_n[i];
                    if(i != 0) {

                        multisin *=  sin2theta_n[i-1];
                    }
                    y2_n[i] = y2 *  multisin * cos2theta_n[i];
                    Eemit[i] = seEpsn_m[seNum-1] * y2_n[i];
                    Eemisum += Eemit[i];


                }

                Eemit[seNum-1] = Eemit[seNum-2] / cos2theta_n[seNum-2] * sin2theta_n[seNum-2];
                Eemisum += Eemit[seNum-1];

                Eemit[0] = Eemisum/seNum;

            } else {// emit seNum particles
                for(int i = 0; i < seNum - 1; i++) {
                    double mu = sePn_m[seNum-1] * (seNum  - i);
                    double nu =  sePn_m[seNum-1];
                    double rand_n = IpplRandom();
                    sin2theta_n[i] = invbetai(rand_n, mu, nu);
                    cos2theta_n[i] = 1.0 - sin2theta_n[i];
                    if(i != 0) {

                        multisin *=  sin2theta_n[i-1];
                    }
                    y2_n[i] = y2 *  multisin * cos2theta_n[i];
                    Eemit[i] = seEpsn_m[seNum-1] * y2_n[i];


                }

                Eemit[seNum-1] = Eemit[seNum-2] / cos2theta_n[seNum-2] * sin2theta_n[seNum-2];

                /*=====================================================================================================*/
            }
        }

    }

    Vector_t z_unit = TriNorm_l;
    Vector_t x_unit = localX - interCoords_l;
    double tmp = sqrt(dot(x_unit,x_unit));
    x_unit /= tmp;
    Vector_t y_unit = cross(z_unit,x_unit);

    size_t lowMark = itsBunch->getLocalNum();
    if (!nEmissionMode) {// emit only 1 particle with larger charge instead of emit seNum secondaries.
        //if (seNum>0) {// we dont delete particles even when seNum==0.
        double gamma_const = Eemit[0] / Physics::m_e/1.0e9 + 1.0;
        double beta_const = sqrt(1.0 - 1.0 / pow(gamma_const, 2.0));
        double P_emitted = gamma_const * beta_const;

        Vector_t P_local;
        Vector_t P_global = (0.0);

        if( Options::ppdebug ) {

            P_global = P_emitted*TriNorm_l;// 1D for parallel plate benchmark

        } else {
            /*==================3D for Furman-Pivi's Model====================*/
            P_local[2] = P_emitted * cos(emiTheta[0]);
            P_local[1] = P_emitted * sin(emiTheta[0]) * sin(emiPhi[0]);
            P_local[0] = P_emitted * sin(emiTheta[0]) * cos(emiPhi[0]);
            P_global = P_local[0]*x_unit + P_local[1]*y_unit +  P_local[2]*z_unit;//Pivi's model
            /*================================================================*/


        }
        itsBunch->create(1);
        itsBunch->R[lowMark] = interCoords_l;

        itsBunch->P[lowMark] = P_global;
        itsBunch->Bin[lowMark] = 0;
        itsBunch->PType[lowMark] = ParticleType::NEWSECONDARY;
        itsBunch->TriID[lowMark] = 0;
        //itsBunch->Q[lowMark] = incQ_l*seNum;// charge of simulation particle will be sum of secondaies
        itsBunch->Q[lowMark] = incQ_l*seyNum;// charge of simulation particle will be multiplied by SEY.
        itsBunch->Ef[lowMark] = Vector_t(0.0);
        itsBunch->Bf[lowMark] = Vector_t(0.0);
        itsBunch->dt[lowMark] = itsBunch->getdT();
        //}



    } else {
        for(size_t i = 0; i < (size_t) seNum; i++) {

            double gamma_const = Eemit[i] / Physics::m_e/1.0e9 + 1.0;
            double beta_const = sqrt(1.0 - 1.0 / pow(gamma_const, 2.0));
            double P_emitted = gamma_const * beta_const;

            Vector_t P_local;
            Vector_t P_global = (0.0);

            if( Options::ppdebug ) {

                P_global = P_emitted*TriNorm_l;// 1D for parallel plate benchmark

            } else {
                /*==================3D for Furman-Pivi's Model====================*/
                P_local[2] = P_emitted * cos(emiTheta[i]);
                P_local[1] = P_emitted * sin(emiTheta[i]) * sin(emiPhi[i]);
                P_local[0] = P_emitted * sin(emiTheta[i]) * cos(emiPhi[i]);
                P_global = P_local[0]*x_unit + P_local[1]*y_unit +  P_local[2]*z_unit;//Pivi's model
                /*================================================================*/


            }
            itsBunch->create(1);
            itsBunch->R[lowMark+i] = interCoords_l;

            itsBunch->P[lowMark+i] = P_global;
            itsBunch->Bin[lowMark+i] = 0;
            itsBunch->PType[lowMark+i] = ParticleType::NEWSECONDARY;
            itsBunch->TriID[lowMark+i] = 0;
            itsBunch->Q[lowMark+i] = incQ_l;
            itsBunch->Ef[lowMark+i] = Vector_t(0.0);
            itsBunch->Bf[lowMark+i] = Vector_t(0.0);
            itsBunch->dt[lowMark+i] = itsBunch->getdT();

        }
    }

    IpplTimings::stopTimer(TPnSec_m);
}


//Vaughan's secondary emission model.
void SecondaryEmissionPhysics::nSec(const double &incEnergy,
                                    const double &cosTheta,
                                    int &seNum, int &seType,
                                    const double &incQ,
                                    const Vector_t &TriNorm,
                                    const Vector_t &inteCoords,
                                    const Vector_t &localX,
                                    PartBunchBase<double, 3> *itsBunch,
                                    double &seyNum, const double &ppVw,
                                    const double &vSeyZero,
                                    const double &vEzero,
                                    const double &vSeyMax,
                                    const double &vEmax,
                                    const double &vKenergy,
                                    const double &vKtheta,
                                    const double &vVThermal,
                                    const bool nEmissionMode)
{


    IpplTimings::startTimer(TPnSec_m);

    std::vector<Vector_t> se_P;
    calcEmiNum(incEnergy, cosTheta, seNum, vSeyZero, vEzero, vSeyMax, vEmax, vKenergy, vKtheta, seyNum);//calculate emitted number & SEY factor
    double Eemit[seNum];
    double emiTheta[seNum];
    double emiPhi[seNum];
    Vector_t interCoords_l = inteCoords;
    Vector_t TriNorm_l = TriNorm;//simpler model
    double vw=ppVw; //1.6*1e-19*1200/9.10938188*1e-31/(2*3.1415926*2.0*1e8)/0.03;//benchmark
    double vt=vVThermal;//7.268929821*1e5 1.5eV//benchmark
    double f_max;//benchmark
    double test_a;//benchmark
    double test_asq;//benchmark
    if( Options::ppdebug ) {
        f_max=vw/vt*exp(-0.5);// velocity a Maxwellian distribution. See Anza et al.,Phys. Plasmas 17, 062110 (2010)

        test_a=vt/vw;
        test_asq=test_a*test_a;
    } else {// Energy Maxwell-Boltzmann distribution f(E)=E/a^2*exp(-E/a), a= kinetic energy w.r.t thermal speed vt. See www.nuc.berkeley.edu/dept/Courses/NE-255/minicourse.pdf p.20
        test_a = Physics::m_e*(1.0/sqrt(1-vt*vt/Physics::c/Physics::c)-1.0)*1.0e9;// m_e GeV change it to eV
        test_asq=test_a*test_a;
        f_max= 1.0/test_a*exp(-1.0);

    }

    double incQ_l = incQ;
    if( Options::ppdebug ) {
        // 1D emission angle along the surface triangle normal

    } else {
        if (!nEmissionMode) {

            double tmp1 = IpplRandom();
            double tmp2 = IpplRandom();
            double seAlpha = 1.0;
            double temp = 1.0 / (1.0 + seAlpha);// pow(cosine(theta),seAlpha) distribution. Here seAlpha=1.0
            emiTheta[0] = acos(pow(tmp1, temp));
            emiPhi[0] = Physics::two_pi * tmp2;


        } else {

            for(int i = 0; i < seNum; i++) {

                double tmp1 = IpplRandom();
                double tmp2 = IpplRandom();
                double seAlpha = 1.0;
                double temp = 1.0 / (1.0 + seAlpha);// pow(cosine(theta),seAlpha) distribution. Here seAlpha=1.0
                emiTheta[i] = acos(pow(tmp1, temp));
                emiPhi[i] = Physics::two_pi * tmp2;

            }

        }
    }

    if(seNum == 0) {
        if (!nEmissionMode) {
            if( Options::ppdebug ) {
                double test_s=1.0;
                double f_x=0;
                double test_x=0;
                while (test_s>f_x) {
                    test_s=IpplRandom();
                    test_s*=f_max;
                    test_x=IpplRandom();
                    test_x*=10*test_a;// truncation range for emission speed(0,10*test_a);
                    f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
                }
                double v_emi=test_x*vw;
                Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/Physics::c/Physics::c)-1)*Physics::m_e*1.0e9;
            } else {
                double test_s=1.0;
                double f_x=0;
                double test_x=0;
                while (test_s>f_x) {// rejection & acceptance method
                    test_s=IpplRandom();
                    test_s*=f_max;
                    test_x=IpplRandom();
                    test_x*=10*test_a;// truncation range for emission energy(0,10*test_a);
                    f_x=test_x/test_asq*exp(-test_x/test_a);// thermal distribution for energy
                }

                Eemit[0]=test_x;


            }

        }
        // else The incident particle will be marked for deletion


    } else {
        seType = 2;
        if (!nEmissionMode) {
            if( Options::ppdebug ) {
                double test_s=1.0;
                double f_x=0;
                double test_x=0;
                while (test_s>f_x) {
                    test_s=IpplRandom();
                    test_s*=f_max;
                    test_x=IpplRandom();
                    test_x*=10*test_a;// truncation range for emission speed(0,10*test_a);
                    f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
                }
                double v_emi=test_x*vw;
                Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/Physics::c/Physics::c)-1)*Physics::m_e*1.0e9;
            } else {
                double test_s=1.0;
                double f_x=0;
                double test_x=0;
                while (test_s>f_x) {
                    test_s=IpplRandom();
                    test_s*=f_max;
                    test_x=IpplRandom();
                    test_x*=10*test_a;// truncation range for emission energy(0,10*test_a);
                    f_x=test_x/test_asq*exp(-test_x/test_a);// thermal distribution for energy
                }

                Eemit[0]=test_x;

            }

        } else {
            if( Options::ppdebug ) {
                /*=======================Maxwellian Distribution====================*/
                // the velocity distribution for Vaughan's model:fu=u*vw*vw/vt/vt*exp(-u*u*vw*vw/vt/vt) valid for benchmarking;

                for(int i = 0; i < seNum; i++) {
                    double test_s=1.0;
                    double f_x=0;
                    double test_x=0;
                    while (test_s>f_x) {
                        test_s=IpplRandom();
                        test_s*=f_max;
                        test_x=IpplRandom();
                        test_x*=10*test_a;//range for normalized emission speed(0,10*test_a);
                        f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
                    }
                    double v_emi=test_x*vw;
                    Eemit[i]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;

                }

                /*---------------------End Of Maxwellian Distribution--------------------*/
            } else {// energy thermal distribution for Vaughan model
                for(int i = 0; i < seNum; i++) {
                    double test_s=1.0;
                    double f_x=0;
                    double test_x=0;
                    while (test_s>f_x) {
                        test_s=IpplRandom();
                        test_s*=f_max;
                        test_x=IpplRandom();
                        test_x*=10*test_a;// truncation range for emission energy(0,10*test_a);
                        f_x=test_x/test_asq*exp(-test_x/test_a);// thermal distribution for energy
                    }

                    Eemit[i]=test_x;
                }
            }
        }
    }

    Vector_t z_unit = TriNorm;
    Vector_t x_unit = localX - interCoords_l;
    double tmp = sqrt(dot(x_unit,x_unit));
    x_unit /= tmp;
    Vector_t y_unit = cross(z_unit,x_unit);

    size_t lowMark = itsBunch->getLocalNum();
    if (!nEmissionMode) {
        double gamma_const = Eemit[0] / Physics::m_e/1.0e9 + 1.0;
        double beta_const = sqrt(1.0 - 1.0 / pow(gamma_const, 2.0));
        double P_emitted = gamma_const * beta_const;


        Vector_t P_local;
        Vector_t P_global = (0.0);
        if( Options::ppdebug ) {

            P_global = P_emitted*TriNorm_l;//1D for parallel plate benchmark

        } else {
            /*==================3D Vaughan' Model==========================================*/
            P_local[2] = P_emitted * cos(emiTheta[0]);
            P_local[1] = P_emitted * sin(emiTheta[0]) * sin(emiPhi[0]);
            P_local[0] = P_emitted * sin(emiTheta[0]) * cos(emiPhi[0]);
            P_global = P_local[0]*x_unit + P_local[1]*y_unit +  P_local[2]*z_unit;// the same cosine distribution as Pivi's model
            /*=============================================================================*/
        }

        itsBunch->create(1);
        itsBunch->R[lowMark] = interCoords_l;
        itsBunch->P[lowMark] = P_global;
        itsBunch->Bin[lowMark] = 0;
        itsBunch->PType[lowMark] = ParticleType::NEWSECONDARY;
        itsBunch->TriID[lowMark] = 0;
        itsBunch->Q[lowMark] = incQ_l*seyNum;
        itsBunch->Ef[lowMark] = Vector_t(0.0);
        itsBunch->Bf[lowMark] = Vector_t(0.0);
        itsBunch->dt[lowMark] = itsBunch->getdT();

    } else {
        for(size_t i = 0; i < (size_t) seNum; i++) {

            double gamma_const = Eemit[i] / Physics::m_e/1.0e9 + 1.0;
            double beta_const = sqrt(1.0 - 1.0 / pow(gamma_const, 2.0));
            double P_emitted = gamma_const * beta_const;


            Vector_t P_local;
            Vector_t P_global = (0.0);
            if( Options::ppdebug ) {

                P_global = P_emitted*TriNorm_l;//1D for parallel plate benchmark

            } else {
                /*==================3D Vaughan' Model==========================================*/
                P_local[2] = P_emitted * cos(emiTheta[i]);
                P_local[1] = P_emitted * sin(emiTheta[i]) * sin(emiPhi[i]);
                P_local[0] = P_emitted * sin(emiTheta[i]) * cos(emiPhi[i]);
                P_global = P_local[0]*x_unit + P_local[1]*y_unit +  P_local[2]*z_unit;//the same cosine distribution as Pivi's model
                /*=============================================================================*/
            }

            itsBunch->create(1);
            itsBunch->R[lowMark+i] = interCoords_l;

            itsBunch->P[lowMark+i] = P_global;
            itsBunch->Bin[lowMark+i] = 0;
            itsBunch->PType[lowMark+i] = ParticleType::NEWSECONDARY;
            itsBunch->TriID[lowMark+i] = 0;
            itsBunch->Q[lowMark+i] = incQ_l;
            itsBunch->Ef[lowMark+i] = Vector_t(0.0);
            itsBunch->Bf[lowMark+i] = Vector_t(0.0);
            itsBunch->dt[lowMark+i] = itsBunch->getdT();


        }
    }

    IpplTimings::stopTimer(TPnSec_m);

}

void SecondaryEmissionPhysics::calcEmiNum(const double incEnergy,
                                          const double cosTheta,
                                          int &seNum,
                                          const double &vSeyZero,
                                          const double &vEzero,
                                          const double &vSeyMax,
                                          const double &vEmax,
                                          const double &vKenergy,
                                          const double &vKtheta,
                                          double &seyNum) {// For Vaughan's model.

    double vSEY = 0;

    if (incEnergy<vEzero) {
        vSEY = vSeyZero;
    } else {

        double theta = acos(cosTheta);
        double delta_max = vSeyMax*(1.0+vKtheta*theta*theta/Physics::two_pi);//here the symbols k are different with reference.
        double E_max = vEmax*(1.0+vKenergy*theta*theta/Physics::two_pi);
        PAssert_GT(E_max - vEzero, 0);
        double v = (incEnergy-vEzero)/(E_max-vEzero);

        if (v<=3.6) {

            if (v<1.0) {
                vSEY = delta_max*pow(v*exp(1.0-v),0.56);
            }else {
                vSEY = delta_max*pow(v*exp(1.0-v),0.25);
            }

        }else {
            vSEY = delta_max*1.125/pow(v,0.35);
        }

    }
    double L = exp(-vSEY);// poisson distribution: Knuth's algorithm.
    int k = 0;
    double p = 1.0;
    do {
        k++;
        double u = IpplRandom();
        p*=u;
    }
    while (p>L);
    seNum = k-1;
    seyNum = vSEY;

}

void SecondaryEmissionPhysics::calcEmiNum(const double /*incEnergy*/, const double /*cosTheta*/, const double *prob, int &seNum) {// For Furman-Pivi's model

    double prob_max = 0.0;
    // Acceptance-rejection methods to generate random number with specified distribution.

    for(int i = 0; i < 11; i++) {

        if(prob[i] > prob_max) {
            prob_max = prob[i];
        }

    }
    double pY  = 1.0;
    double pX = 0.0;
    while(pY > pX) {

        double rand1 = IpplRandom();
        //double rand1 = (*rand)(rng);
        seNum = (int)(rand1 * 11.0);//fix me
        pX = prob[seNum];
        double rand2 = IpplRandom();
        //double rand2 = (*rand)(rng);
        pY = prob_max * rand2;
    }

}


double SecondaryEmissionPhysics::calcDeltats(const double incEnergy, const double cosTheta) {

    double seypeak = seYPeakTS_m * (1 + seTOneTS_m * (1.0 - pow(cosTheta, seTTwoTS_m))); //formula III.E (48a)
    double seepeak = seEPeakTS_m * (1 + seTThreeTS_m * (1.0 - pow(cosTheta, seTFourTS_m))); //formula III.E (48b)
    double tmpx = incEnergy / seepeak;
    double tmpD = seSTS_m * tmpx / (seSTS_m - 1 + pow(tmpx, seSTS_m)); //formula III.D (32)
    double ret = seypeak * tmpD; //formula III.D (31)
    return ret;

}

double SecondaryEmissionPhysics::calcDeltar(const double incEnergy, const double cosTheta) {

    double tmp = pow(incEnergy / seERed_m, seR_m);
    double ret = sePRed_m * (1.0 - exp(-1 * tmp)); //formula III.D (28)
    ret = ret * (1.0 + seROne_m * (1.0 - pow(cosTheta, seRTwo_m))); //formula III.E (47b)
    return ret;

}


double SecondaryEmissionPhysics::calcDeltae(const double incEnergy, const double cosTheta) {

    double tmp = pow(std::abs(incEnergy - seEScatPeak_m) / seW_m, seP_m) / seP_m;
    double ret = sePScat_m + (sePScatPeak_m - sePScat_m) * exp(-1 * tmp); //formula III.D (25)
    ret = ret * (1.0 + seEOne_m * (1.0 - pow(cosTheta, seETwo_m))); //formula III.E (47a)
    return ret;

}


double SecondaryEmissionPhysics::calcProb(const double incEnergy, const double cosTheta, double *prob) {

    deltae_m = calcDeltae(incEnergy, cosTheta);
    deltar_m = calcDeltar(incEnergy, cosTheta);
    deltats_m = calcDeltats(incEnergy, cosTheta);

    double tmp = 1.0 - deltae_m - deltar_m;
    double p = deltats_m / tmp / 10.0;
    double q = 1.0 - p;
    double b[11];
    b[0]  = 1.0;
    b[1]  = 10.0;
    b[2]  = 45.0;
    b[3]  = 120.0;
    b[4]  = 210.0;
    b[5]  = 252.0;
    b[6]  = 210.0;
    b[7]  = 120.0;
    b[8]  = 45.0;
    b[9]  = 10.0;
    b[10] = 1.0;
    for(int i = 0; i < 11; i++) {
        prob[i] = tmp * b[i] * pow(p, i) * pow(q, (10 - i));
    }
    prob[1] = prob[1] + deltae_m + deltar_m;

    /*==============================================*/

    return (deltae_m+deltar_m+deltats_m);
    //cout << "sum prob: " << sum << endl;
    /*==============================================*/
}

void SecondaryEmissionPhysics::setSeMaterial(int material_num) {

    if(material_num == 0) {
        seAlpha_m = 1.0;
        sePScat_m = 0.02;
        sePScatPeak_m = 0.496;
        seEScatPeak_m = 0;
        seW_m = 60.86;
        seP_m = 1.0;
        seDelta_m = 2.0;
        seEOne_m = 0.26;
        seETwo_m = 2;

        sePRed_m = 0.2;
        seERed_m = 0.041;
        seR_m = 0.104;
        seQ_m = 0.5;
        seROne_m = 0.26;
        seRTwo_m = 2;

        seYPeakTS_m = 1.8848;
        seEPeakTS_m = 276.8;
        seSTS_m = 1.54;
        seTOneTS_m = 0.66;
        seTTwoTS_m = 0.8;
        seTThreeTS_m = 0.7;
        seTFourTS_m = 1.0;
        seEPeakTot_m = 271;
        seYPeakTot_m = 2.1;
        sePn_m[0] = 2.5;
        sePn_m[1] = 3.3;
        sePn_m[2] = 2.5;
        sePn_m[3] = 2.5;
        sePn_m[4] = 2.8;
        sePn_m[5] = 1.3;
        sePn_m[6] = 1.5;
        sePn_m[7] = 1.5;
        sePn_m[8] = 1.5;
        sePn_m[9] = 1.5;
        seEpsn_m[0] = 1.5;
        seEpsn_m[1] = 1.75;
        seEpsn_m[2] = 1.0;
        seEpsn_m[3] = 3.75;
        seEpsn_m[4] = 8.5;
        seEpsn_m[5] = 11.5;
        seEpsn_m[6] = 2.5;
        seEpsn_m[7] = 3.0;
        seEpsn_m[8] = 2.5;
        seEpsn_m[9] = 3.0;
    }
    if(material_num == 1) {
        seAlpha_m = 1.0;
        sePScat_m = 0.07;
        sePScatPeak_m = 0.5;
        seEScatPeak_m = 0;
        seW_m = 100;
        seP_m = 0.9;
        seDelta_m = 1.9;
        seEOne_m = 0.26;
        seETwo_m = 2;

        sePRed_m = 0.74;
        seERed_m = 40;
        seR_m = 1;
        seQ_m = 0.4;
        seROne_m = 0.26;
        seRTwo_m = 2;

        seYPeakTS_m = 1.22;
        seEPeakTS_m = 310;
        seSTS_m = 1.813;
        seTOneTS_m = 0.66;
        seTTwoTS_m = 0.8;
        seTThreeTS_m = 0.7;
        seTFourTS_m = 1.0;
        seEPeakTot_m = 292;
        seYPeakTot_m = 2.05;

        sePn_m[0] = 1.6;
        sePn_m[1] = 2.0;
        sePn_m[2] = 1.8;
        sePn_m[3] = 4.7;
        sePn_m[4] = 1.8;
        sePn_m[5] = 2.4;
        sePn_m[6] = 1.8;
        sePn_m[7] = 1.8;
        sePn_m[8] = 2.3;
        sePn_m[9] = 1.8;
        seEpsn_m[0] = 3.9;
        seEpsn_m[1] = 6.2;
        seEpsn_m[2] = 13.0;
        seEpsn_m[3] = 8.8;
        seEpsn_m[4] = 6.25;
        seEpsn_m[5] = 2.25;
        seEpsn_m[6] = 9.2;
        seEpsn_m[7] = 5.3;
        seEpsn_m[8] = 17.8;
        seEpsn_m[9] = 10;
    }
}



/*==========================================
 *http://www.taygeta.com/random/gaussian.html
 *return a gaussian distributed random number
 *==========================================*/

double SecondaryEmissionPhysics::gaussRand() {

    double x1;
    double x2;
    double w;
    do {
        x1 = 2.0 * IpplRandom() - 1.0;
        //x1 = 2.0 * (*rand)(rng) - 1.0;
        x2 = 2.0 * IpplRandom() - 1.0;
        //x2 = 2.0 * (*rand)(rng) - 1.0;
        w = x1 * x1 + x2 * x2;
    } while(w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    double ret = x1 * w;
    return ret;

}

double SecondaryEmissionPhysics::gammp(const double a, const double x) {
    // Returns the incomplete gamma function P .a; x/.
    static const int ASWITCH = 100; //When to switch to quadrature method.
    if(x < 0.0 || a <= 0.0)
        throw("bad args in gammp");
    if(x == 0.0)
        return 0.0;
    else if((int)a >= ASWITCH)
        return gammpapprox(a, x, 1); //Quadrature.
    else if(x < a + 1.0)
        return gser(a, x); //Use the series representation.
    else
        return 1.0 - gcf(a, x); //Use the continued fraction representation.
}

double SecondaryEmissionPhysics::gser(const double a, const double x) {
    // Returns the incomplete gamma function P .a; x/ evaluated by its series representation.
    // Also sets ln .a/ as gln. User should not call directly.
    double sum, del, ap;
    double gln = gammln(a);
    ap = a;
    del = sum = 1.0 / a;
    for(;;) {
        ++ap;
        del *= x / ap;
        sum += del;
        if(std::abs(del) < std::abs(sum)*myeps::EPS) {
            return sum * exp(-x + a * log(x) - gln);
        }
    }
}

double SecondaryEmissionPhysics::gcf(const double a, const double x) {
    // Returns the incomplete gamma function Q.a; x/ evaluated by its continued fraction representation. Also sets ln .a/ as gln. User should not call directly.
    int i;
    double an, b, c, d, del, h;
    double gln = gammln(a);
    b = x + 1.0 - a; // Set up for evaluating continued fraction
    c = 1.0 / myeps::FPMIN; // by modified Lentz¡¯s method (5.2)
    d = 1.0 / b; // with b0 D 0.
    h = d;
    for(i = 1;; i++) {
        //Iterate to convergence.
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if(std::abs(d) < myeps::FPMIN) d = myeps::FPMIN;
        c = b + an / c;
        if(std::abs(c) < myeps::FPMIN) c = myeps::FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if(std::abs(del - 1.0) <= myeps::EPS) break;
    }
    return exp(-x + a * log(x) - gln) * h; //Put factors in front.
}

double SecondaryEmissionPhysics::gammpapprox(double a, double x, int psig) {
    // Incomplete gamma by quadrature. Returns P .a; x/ or Q.a; x/, when psig is 1 or 0,respectively. User should not call directly.

    const double y[18] = {0.0021695375159141994,
                          0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
                          0.082502225484340941, 0.12007019910960293, 0.16415283300752470,
                          0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
                          0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
                          0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
                          0.87126389619061517, 0.95698180152629142
    };
    const double w[18] = {0.0055657196642445571,
                          0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
                          0.034213810770299537, 0.040875750923643261, 0.047235083490265582,
                          0.053244713977759692, 0.058860144245324798, 0.064039797355015485,
                          0.068745323835736408, 0.072941885005653087, 0.076598410645870640,
                          0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
                          0.085346685739338721, 0.085983275670394821
    };
    int j;
    double xu, t, sum, ans;
    double a1 = a - 1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
    int ngau = 18;
    double gln = gammln(a);// Set how far to integrate into the tail:
    if(x > a1) xu = max_Gamma(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
    else xu = max_Gamma(0., min_Gamma(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));
    sum = 0;
    for(j = 0; j < ngau; j++) {
        //Gauss-Legendre.
        t = x + (xu - x) * y[j];
        sum += w[j] * exp(-(t - a1) + a1 * (log(t) - lna1));
    }
    ans = sum * (xu - x) * exp(a1 * (lna1 - 1.) - gln);
    return (psig ? (ans > 0.0 ? 1.0 - ans : -ans) : (ans >= 0.0 ? ans : 1.0 + ans));
}

double SecondaryEmissionPhysics::gammln(const double xx) {
    // Returns the value ln%GÅ%@.xx/ for xx > 0.
    int j;
    double x, tmp, y, ser;
    static const double cof[14] = {57.1562356658629235, -59.5979603554754912,
                                   14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
                                   .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,
                                   -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
                                   .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5
    };
    if(xx <= 0) throw("bad arg in gammln");
    y = x = xx;
    tmp = x + 5.24218750000000000; // Rational 671/128.
    tmp = (x + 0.5) * log(tmp) - tmp;
    ser = 0.999999999999997092;
    for(j = 0; j < 14; j++) ser += cof[j] / ++y;
    return tmp + log(2.5066282746310005 * ser / x);
}

double SecondaryEmissionPhysics::invgammp(double p, double a) {
    //Returns x such that P .a; x/ D p for an argument p between 0 and 1.
    int j;
    double x, err, t, u, pp, lna1 = 0.0, afac = 0.0, a1 = a - 1;
    const double EPS = 1.e-8; //Accuracy is the square of EPS.
    double gln = gammln(a);
    if(a <= 0.) throw("a must be pos in invgammap");
    if(p >= 1.) return max_Gamma(100., a + 100.*sqrt(a));
    if(p <= 0.) return 0.0;
    if(a > 1.) {
        //Initial guess based on reference [1].
        lna1 = log(a1);
        afac = exp(a1 * (lna1 - 1.) - gln);
        pp = (p < 0.5) ? p : 1. - p;
        t = sqrt(-2.*log(pp));
        x = (2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t;
        if(p < 0.5) x = -x;
        x = max_Gamma(1.e-3, a * pow(1. - 1. / (9.*a) - x / (3.*sqrt(a)), 3));
    } else {
        //Initial guess based on equations (6.2.8) and (6.2.9).
        t = 1.0 - a * (0.253 + a * 0.12);

        if(p < t) x = pow(p / t, 1. / a);
        else x = 1. - log(1. - (p - t) / (1. - t));
    }
    for(j = 0; j < 12; j++) {
        if(x <= 0.0) return 0.0;
        //x too small to compute accurately.
        err = gammp(a, x) - p;
        if(a > 1.) t = afac * exp(-(x - a1) + a1 * (log(x) - lna1));
        else t = exp(-x + a1 * log(x) - gln);
        u = err / t;
        x -= (t = u / (1. - 0.5 * min_Gamma(1., u * ((a - 1.) / x - 1))));
        //Halley¡¯s method.
        if(x <= 0.) x = 0.5 * (x + t);
        //Halve old value if x tries to go negative.
        if(std::abs(t) < EPS * x) break;
    }
    return x;
}

double SecondaryEmissionPhysics::betai(const double x, const double a, const double b) {
    //cout<<"betai called"<<endl;
    //Returns incomplete beta function Ix.a; b/ for positive a and b, and x between 0 and 1.
    static const int SWITCH = 3000;
    double bt;
    /*=========================debug code=====================
      if (a==b) {
      cout<< " x in betai "<<x<<endl;
      }
      =========================debug code=====================*/
    if(a <= 0.0 || b <= 0.0) {
        //cout<<"betai 1"<<endl;
        throw("Bad a or b in routine betai");
    }

    if(x < 0.0 || x > 1.0) {
        //cout<<"betai 2"<<endl;
        throw("Bad x in routine betai");

    }
    if(x == 0.0 || x == 1.0) {
        //cout<<"betai 3"<<endl;
        return x;
    }
    if(a > SWITCH && b > SWITCH) {
        //cout<<"betai 4"<<endl;
        return betaiapprox(a, b, x);
    }
    bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
    if(x < (a + 1.0) / (a + b + 2.0)) {
        //cout<<"betai 5"<<endl;
        return bt * betacf(a, b, x) / a;
    } else {
        //cout<<"betai 6"<<endl;
        /*=========================debug code=====================
          if (a==b) {
          cout<<" betacf(b, a, 1.0 - x) = "<< betacf(b, a, 1.0 - x)<<" a = "<<a<<" b = "<<b<<endl;
          }
          =========================================================*/
        return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
    }
}

double SecondaryEmissionPhysics::betacf(const double a, const double b, const double x) {
    //Evaluates continued fraction for incomplete beta function by modified Lentzâs method (5.2). User should not call directly.
    int m, m2;
    double aa, c, d, del, h, qab, qam, qap;
    qab = a + b; //These qab will be used in factors that occur in the coefficients (6.4.6).
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0; //First step of Lentzâs method.

    d = 1.0 - qab * x / qap;
    /*========debug====================
      if (a == b) {
      cout<<"d = "<<d<<" qap = "<<qap<<" qab = "<<qab<<" x = "<<x<<endl;
      }
      ========debug====================*/
    if(std::abs(d) < myeps::FPMIN) d = myeps::FPMIN;
    d = 1.0 / d;
    h = d;
    // cout<<"h = "<<h<<" myeps::FPMIN = "<<myeps::FPMIN<<endl;
    for(m = 1; m < 10000; m++) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d; //One step (the even one) of the recurrence.
        if(std::abs(d) < myeps::FPMIN)
            d = myeps::FPMIN;
        c = 1.0 + aa / c;
        if(std::abs(c) < myeps::FPMIN)
            c = myeps::FPMIN;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d; //Next step of the recurrence (the odd one).
        if(std::abs(d) < myeps::FPMIN)
            d = myeps::FPMIN;
        c = 1.0 + aa / c;
        if(std::abs(c) < myeps::FPMIN)
            c = myeps::FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if(std::abs(del - 1.0) <= myeps::EPS)
            break;
    }
    return h;
}

double SecondaryEmissionPhysics::betaiapprox(double a, double b, double x) {
    //Incomplete beta by quadrature. Returns Ix.a; b/. User should not call directly.
    const double y[18] = {0.0021695375159141994,
                          0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
                          0.082502225484340941, 0.12007019910960293, 0.16415283300752470,
                          0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
                          0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
                          0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
                          0.87126389619061517, 0.95698180152629142
    };
    const double w[18] = {0.0055657196642445571,
                          0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
                          0.034213810770299537, 0.040875750923643261, 0.047235083490265582,
                          0.053244713977759692, 0.058860144245324798, 0.064039797355015485,
                          0.068745323835736408, 0.072941885005653087, 0.076598410645870640,
                          0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
                          0.085346685739338721, 0.085983275670394821
    };
    int j;
    double xu, t, sum, ans;
    double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
    double lnmu = log(mu), lnmuc = log(1. - mu);
    t = sqrt(a * b / (sqrt(a + b) * (a + b + 1.0)));
    if(x > a / (a + b)) { //Set how far to integrate into the tail:
        if(x >= 1.0)
            return 1.0;
        xu = min_Gamma(1., max_Gamma(mu + 10.*t, x + 5.0 * t));
    } else {

        if(x <= 0.0)
            return 0.0;
        xu = max_Gamma(0., min_Gamma(mu - 10.*t, x - 5.0 * t));
    }
    sum = 0;
    for(j = 0; j < 18; j++) { //Gauss-Legendre.
        t = x + (xu - x) * y[j];
        sum += w[j] * exp(a1 * (log(t) - lnmu) + b1 * (log(1 - t) - lnmuc));
    }
    ans = sum * (xu - x) * exp(a1 * lnmu - gammln(a) + b1 * lnmuc - gammln(b) + gammln(a + b));
    return ans > 0.0 ? 1.0 - ans : -ans;
}

double SecondaryEmissionPhysics::invbetai(double p, double a, double b) {
    // Inverse of incomplete beta function. Returns x such that Ix.a; b/ D p for argument p between 0 and 1.
    const double EPS = 1.e-8;
    double pp, t, u, err, x, al, h, w, afac, a1 = a - 1., b1 = b - 1.;
    int j;

    if(p <= 0.)
        return 0.;
    else if(p >= 1.)
        return 1.;
    else if(a >= 1. && b >= 1.) { // Set initial guess. See text.
        pp = (p < 0.5) ? p : 1. - p;
        t = sqrt(-2.*log(pp));
        x = (2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t;
        /*========debug====================
          if (a==b) {

          cout<<"x = "<<x<<" p = "<<p<<" w = "<<w<<" t = "<<t<<" al = "<<al<<" h = "<<h<<endl;

          }
          ========debug====================*/
        // if(p < 0.5)//origin code from numerical ricipes.
        if(x < 0.0)//fixed bug.
            x = -x;
        al = (sqrt(x) - 3.) / 6.;
        h = 2. / (1. / (2.*a - 1.) + 1. / (2.*b - 1.));
        w = (x * sqrt(al + h) / h) - (1. / (2.*b - 1) - 1. / (2.*a - 1.)) * (al + 5. / 6. - 2. / (3.*h));
        x = a / (a + b * exp(2.*w));

    } else {
        double lna = log(a / (a + b)), lnb = log(b / (a + b));
        t = exp(a * lna) / a;
        u = exp(b * lnb) / b;
        w = t + u;
        if(p < t / w)
            x = pow(a * w * p, 1. / a);
        else
            x = 1. - pow(b * w * (1. - p), 1. / b);
    }
    afac = -gammln(a) - gammln(b) + gammln(a + b);
    for(j = 0; j < 10; j++) {
        if(x == 0. || x == 1.) // a or b too small for accurate calculation.
            return x;
        err = betai(x, a, b) - p;
        /*=========================debug code=====================
          if (a==b) {
          cout<<"p = "<<p<<" err = "<<err<<" t = "<<t<<" u = "<<u<<" x = "<<x<<endl;
          }
          =========================================================*/
        t = exp(a1 * log(x) + b1 * log(1. - x) + afac);
        u = err / t; //Halley:
        x -= (t = u / (1. - 0.5 * min_Gamma(1., u * (a1 / x - b1 / (1. - x)))));
        /*=========================debug code=====================
          if (a==b) {
          cout<<"err = "<<err<<" t = "<<t<<" u = "<<u<<" x = "<<x<<endl;
          }
          =========================================================*/
        if(x <= 0.)
            x = 0.5 * (x + t); // Bisect if x tries to go neg or > 1.
        if(x >= 1.)
            x = 0.5 * (x + t + 1.);
        if(std::abs(t) < EPS * x && j > 0)
            break;
    }
    return x;
}