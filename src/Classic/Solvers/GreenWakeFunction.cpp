#include "Solvers/GreenWakeFunction.hh"
#include "Algorithms/PartBunchBase.h"
#include "Utilities/GeneralClassicException.h"

#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <iostream>  // Needed for stream I/O
#include <iomanip>   // Needed for I/O manipulators
#include "gsl/gsl_fft_real.h"
#include "gsl/gsl_fft_halfcomplex.h"

using namespace std;

//IFF: TEST
//#define ENABLE_WAKE_DEBUG
//#define ENABLE_WAKE_DUMP
//#define ENABLE_WAKE_TESTS_FFT_OUT


/**
 * @brief   just a Testfunction!  Calculate the energy of the Wakefunction with the lambda
 *
 *
 * @todo        In this code one can only apply either the longitudinal wakefield or the transversal wakefield. One should implement that both wakefields can be applied to the particle beam
 * @todo        NBins must be set equal to MT of the fieldsolver. This should be changed. (the length of lineDensity_m must be NBins and not to MT)
 *
 * @param[in]   ref
 * @param[in]   NBIN number of Bins
 * @param[in]   Z0 impedance of the tube
 * @param[in]   radius radius of the tube
 * @param[in]   sigma material constant
 * @param[in]   acMode 1 for AC and 2 for DC
 * @param[in]   tau material constant
 * @param[in]   direction 0 for transversal and 1 for Longitudinal
 * @param[in]   constLength  true if the length of the particle bunch is considered as constant
 * @param[in]   fname read wake from file
 */
GreenWakeFunction::GreenWakeFunction(const std::string &name,
                                     ElementBase *element,
                                     vector<Filter *> filters,
                                     int NBIN,
                                     double Z0,
                                     double radius,
                                     double sigma,
                                     int acMode,
                                     double tau,
                                     int direction,
                                     bool constLength,
                                     std::string fname):
    WakeFunction(name, element, NBIN),
    lineDensity_m(),
    //~ FftWField_m(0),
    NBin_m(NBIN),
    Z0_m(Z0),
    radius_m(radius),
    sigma_m(sigma),
    acMode_m(acMode),
    tau_m(tau),
    direction_m(direction),
    constLength_m(constLength),
    filename_m(fname),
    filters_m(filters.begin(), filters.end()) {
#ifdef ENABLE_WAKE_DEBUG
    *gmsg << "* ************* W A K E ************************************************************ " << endl;
    *gmsg << "* Entered GreenWakeFunction::GreenWakeFunction " << '\n';
    *gmsg << "* ********************************************************************************** " << endl;
#endif
}

GreenWakeFunction::~GreenWakeFunction() {
    //~ if(FftWField_m != 0) {
        //~ delete[] FftWField_m;
    //~ }
}


/**
 * @brief   given a vector of length N, distribute the indexes among the available processors
 *
 *
 * @todo        make this function general available
 *
 *
 * @param[in]   length of vector
 * @param[out]  first: lowIndex, second: hiIndex
 */
pair<int, int> GreenWakeFunction::distrIndices(int vectLen) {

    pair<int, int> dist;

    //IFF: properly distribute remainder
    int rem = vectLen - (vectLen / Ippl::getNodes()) * Ippl::getNodes();
    int tmp = (rem > Ippl::myNode()) ? 1 : 0;
    int locBunchRange = vectLen / Ippl::getNodes() + tmp;

    dist.first = locBunchRange * Ippl::myNode() + (1 - tmp) * rem;
    dist.second = dist.first + locBunchRange - 1;

    return dist;
}

void GreenWakeFunction::apply(PartBunchBase<double, 3> *bunch) {

    Vector_t rmin, rmax;
    double charge = bunch->getChargePerParticle();
    // CKR: was bunch->Q[1] changed it;
    //FIXME: why 1? bunch,getTotalCharge()
    // or bunch->getChargePerParticle()?
    double K = 0; // constant to normalize the lineDensity_m to 1
    double spacing, mindist;
    std::vector<double> OutEnergy(NBin_m);

    bunch->calcBeamParameters();
    bunch->get_bounds(rmin, rmax);
    //FIXME IFF: do we have unitless r's here? is that what we want?

    mindist = rmin(2);
    switch(direction_m) {
        case LONGITUDINAL:
            spacing = abs(rmax(2) - rmin(2));
            break; //FIXME: Kann mann das Spacing immer ändern?
        case TRANSVERSAL:
            spacing = rmax(0) * rmax(0) + rmax(1) * rmax(1);
            break;
        default:
            throw GeneralClassicException("GreenWakeFunction", "invalid direction specified");
    }
    assert(NBin_m > 0);
    spacing /= (NBin_m - 1); //FIXME: why -1? CKR: because grid spacings = grid points - 1

    // Calculate the Wakefield if needed
    if(FftWField_m.empty()) {
        FftWField_m.resize(2*NBin_m-1);
        if(filename_m != "") {
            setWakeFromFile(NBin_m, spacing);
        } else {
            CalcWakeFFT(spacing);
        }
    } else if(!constLength_m) {
        CalcWakeFFT(spacing);
    }

    // Calculate the line density of the particle bunch
    std::pair<double, double> meshInfo;
    bunch->calcLineDensity(nBins_m, lineDensity_m, meshInfo);

#ifdef ENABLE_WAKE_DEBUG
    *gmsg << "* ************* W A K E ************************************************************ " << endl;
    *gmsg << "* GreenWakeFunction::apply  lineDensity_m.size() = " << lineDensity_m.size() << endl;
    *gmsg << "* ********************************************************************************** " << endl;
#endif

    // smooth the line density of the particle bunch
    for(vector<Filter *>::const_iterator fit = filters_m.begin(); fit != filters_m.end(); ++fit) {
        (*fit)->apply(lineDensity_m);
    }

    for(unsigned int i = 0; i < lineDensity_m.size(); i++) {
        K += lineDensity_m[i];
    }
    K = 1 / K;

    // compute the kick due to the wakefield
    compEnergy(K, charge, lineDensity_m, OutEnergy.data());

    // Add the right OutEnergy[i] to all the particles
    //FIXME: can we specify LONG AND TRANS?
    switch(direction_m) {
        case LONGITUDINAL:
            for(unsigned int i = 0; i < bunch->getLocalNum(); i++) {

                //FIXME: Stimmt das????????? (von den einheiten)
                // calculate bin containing particle
                int idx = (int)(floor((bunch->R[i](2) - mindist) / spacing));
                //IFF: should be ok
                if(idx == NBin_m) idx--;
                assert(idx >= 0 && idx < NBin_m);
                double dE = OutEnergy[idx];
                bunch->Ef[i](2) += dE;

            }
            break;

        case TRANSVERSAL:
            for(unsigned int i = 0; i < bunch->getLocalNum(); i++) {

                // calculate bin containing particle
                int idx = (int)(floor((bunch->R[i](2) - mindist) / spacing));
                //IFF: should be ok
                if(idx == NBin_m) idx--;
                assert(idx >= 0 && idx < NBin_m);
                double dE = OutEnergy[idx];

                // ACHTUNG spacing auch in transversal richtung
                double dist = sqrt(bunch->R[i](0) * bunch->R[i](0) + bunch->R[i](1) * bunch->R[i](1));
                assert(dist > 0);
                bunch->Ef[i](0) += dE * bunch->R[i](0) / dist;
                bunch->Ef[i](1) += dE * bunch->R[i](1) / dist;

            }
            break;

        default:
            throw GeneralClassicException("GreenWakeFunction", "invalid direction specified");
    }

#ifdef ENABLE_WAKE_DUMP
    ofstream  f2("OutEnergy.dat");
    f2 << "# Energy of the Wake calculated in Opal\n"
       << "# Z0 = " << Z0_m << "\n"
       << "# radius = " << radius << "\n"
       << "# sigma = " << sigma << "\n"
       << "# c = " << c << "\n"
       << "# acMode = " << acMode << "\n"
       << "# tau = " << tau << "\n"
       << "# direction = " << direction << "\n"
       << "# spacing = " << spacing << "\n"
       << "# Lbunch = " << NBin_m << "\n";
    for(int i = 0; i < NBin_m; i++) {
        f2 << i + 1 << " " << OutEnergy[i] << "\n";
    }
    f2.flush();
    f2.close();
#endif
}

/**
 * @brief   just a Testfunction!  Calculate the energy of the Wakefunction with the lambda
 *
 *
 * @param[in]   K a constant
 * @param[in]   charge a constant
 * @param[in]   lambda the distribution of the Particles
 * @param[out]  OutEnergy this is the Output
 */
void GreenWakeFunction::compEnergy(const double K,
                                   const double charge,
                                   const double *lambda,
                                   double *OutEnergy) {
    int N = 2 * NBin_m - 1;
    // Allocate Space for the zero padded lambda and its Fourier Transformed
    std::vector<double> pLambda(N);

    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(N);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(N);

    // fill the arrays with data
    for(int i = 0; i < NBin_m; i ++) {
        pLambda[N-i-1] = 0.0;
        pLambda[i] = lambda[i];
    }

    //FFT of the lambda
    gsl_fft_real_transform(pLambda.data(), 1, N, real, work);
    gsl_fft_real_wavetable_free(real);

    // convolution -> multiplication in Fourier space
    pLambda[0] *= FftWField_m[0];
    for(int i = 1; i < N; i += 2) {
        double temp = pLambda[i];
        pLambda[i] = FftWField_m[i] * pLambda[i] - FftWField_m[i+1] * pLambda[i+1];
        pLambda[i+1] = FftWField_m[i] * pLambda[i+1] + FftWField_m[i+1] * temp;
    }

    // inverse transform to get c, the convolution of a and b;
    hc = gsl_fft_halfcomplex_wavetable_alloc(N);

    gsl_fft_halfcomplex_inverse(pLambda.data(), 1, N, hc, work);

    // Write the result to the output:
    for(int i = 0; i < NBin_m; i ++) {
        OutEnergy[i] = -charge * K * pLambda[i] / (2.0 * NBin_m) * N; // CKR: I doubt that the multiplication with N is correct,
        //      put it here to get the same result as with FFTW
        //      My suspicion: S. Pauli has forgotten that if you
        //      do an fft followed by an inverse fft you'll get
        //      N times your original data

    }


    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);
}

/**
 * @brief   Calculate the energy of the Wakefunction with the lambda
 *
 *
 * @param[in]   K a constant
 * @param[in]   charge a constant
 * @param[in]   lambda the distribution of the Particles
 * @param[out]  OutEnergy this is the Output
 *
 */
void GreenWakeFunction::compEnergy(const double K,
                                   const double charge,
                                   vector<double> lambda,
                                   double *OutEnergy) {
    int N = 2 * NBin_m - 1;
    // Allocate Space for the zero padded lambda and its Fourier Transformed
    std::vector<double> pLambda(N);

    gsl_fft_halfcomplex_wavetable *hc;
    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(N);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(N);

    // fill the arrays with data
    for(int i = 0; i < NBin_m; i ++) {
        pLambda[N-i-1] = 0.0;
        pLambda[i] = lambda[i];
    }

    //FFT of the lambda
    gsl_fft_real_transform(pLambda.data(), 1, N, real, work);
    gsl_fft_real_wavetable_free(real);


    // Convolution -> just a multiplication in Fourier space
    pLambda[0] *= FftWField_m[0];
    for(int i = 1; i < N; i += 2) {
        double temp = pLambda[i];
        pLambda[i] = FftWField_m[i] * pLambda[i] - FftWField_m[i+1] * pLambda[i+1];
        pLambda[i+1] = FftWField_m[i] * pLambda[i+1] + FftWField_m[i+1] * temp;
    }

    // IFFT
    hc = gsl_fft_halfcomplex_wavetable_alloc(N);

    gsl_fft_halfcomplex_inverse(pLambda.data(), 1, N, hc, work);

    // Write the result to the output:
    for(int i = 0; i < NBin_m; i ++) {
        OutEnergy[i] = -charge * K * pLambda[i] / (2.0 * NBin_m) * N; // CKR: I doubt that the multiplication with N is correct,
        //      put it here to get the same result as with FFTW
        //      My suspicion: S. Pauli has forgotten that if you
        //      do an fft followed by an inverse fft you'll get
        //      N times your original data

    }


    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(work);
}

/**
 * @brief   Calculate the FFT of the Wakefunction
 *
 * @param[in]   spacing distance between 2 slice in the line distribution
 *
 */
void GreenWakeFunction::CalcWakeFFT(double spacing) {
    // Set integration properties
    double  a = 1, b = 1000000;
    unsigned int N = 1000000;
    int M = 2 * NBin_m - 1;

    gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(M);
    gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(M);

    const pair<int, int> myDist = distrIndices(NBin_m);
    const int lowIndex = myDist.first;
    const int hiIndex  = myDist.second;

    for(int i = 0; i < M; i ++) {
        FftWField_m[i] = 0.0;
    }

    /**
      Calculate the Wakefield on all processors
      */

    //     if (Ippl::myNode() != Ippl::getNodes()-1) {
    for(int i = lowIndex; i <= hiIndex; i ++) {
        Wake w(i * spacing, Z0_m, radius_m, sigma_m, acMode_m, tau_m, direction_m);
        FftWField_m[i] = simpson(w, a, b, N);
    }
    //     } else {
    //         //IFF: changed  to <= with new distr
    //         for (int i = lowIndex; i <= hiIndex; i ++) {
    //             Wake w(i*spacing, Z0_m, radius_m, sigma_m, acMode_m, tau_m, direction_m);
    //             FftWField[i] = simpson(w,a,b,N);
    //         }
    //     }

    /**
      Reduce the results
      */
    reduce(&(FftWField_m[0]), &(FftWField_m[0]) + NBin_m, &(FftWField_m[0]), OpAddAssign());

#ifdef ENABLE_WAKE_TESTS_FFT_OUT
    std::vector<double> wf(2*NBin_m-1);
    for(int i = 0; i < 2 * NBin_m - 1; ++ i) {
        wf[i] = FftWField_m[i];
    }
#endif
    // calculate the FFT of the Wakefield
    gsl_fft_real_transform(FftWField_m.data(), 1, M, real, work);


#ifdef ENABLE_WAKE_TESTS_FFT_OUT
    ofstream  f2("FFTwake.dat");
    f2 << "# FFT of the Wake calculated in Opal" << "\n"
       << "# Z0 = " << Z0_m << "\n"
       << "# radius = " << radius_m << "\n"
       << "# sigma = " << sigma_m << "\n"
       << "# mode = " << acMode_m << "\n"
       << "# tau = " << tau_m << "\n"
       << "# direction = " << direction_m << "\n"
       << "# spacing = " << spacing << "\n"
       << "# Lbunch = " << NBin_m << "\n";

    f2 << "0\t" << FftWField_m[0] << "\t0.0\t" << wf[0] << "\n";
    for(int i = 1; i < M; i += 2) {
        f2 << (i + 1) / 2 << "\t"
           << FftWField_m[i] << "\t"
           << FftWField_m[i + 1] << "\t"
           << wf[(i+1)/2] << "\n";
    }


    f2.flush();
    f2.close();
#endif
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
}

/**
 * @brief   reads in the wakefield from file
 */
void GreenWakeFunction::setWakeFromFile(int NBin_m, double spacing) {
    Inform msg("readSDDS ");
    std::string name;
    char temp[256];
    int Np;
    double dummy;
    gsl_fft_real_wavetable *real;
    gsl_fft_real_workspace *work;
    std::ifstream fs;

    fs.open(filename_m.c_str());

    if(fs.fail()) {
        throw GeneralClassicException("GreenWakeFunction::setWakeFromFile",
                            "Open file operation failed, please check if \""
                            + filename_m +  "\" really exists.");
    }

    fs >> name;
    msg << " SSDS1 read = " << name << endl;
    if(name.compare("SDDS1") != 0) {
        throw GeneralClassicException("GreenWakeFunction::setWakeFromFile",
                            " No SDDS1 File. A SDDS1 file should start with a SDDS1 String. Check file \""
                            + filename_m +  "\" ");
    }

    for(int i = 0; i < 6; i++) {
        fs.getline(temp, 256);
        msg << "line " << i << " :   " << temp << endl;
    }

    fs >> Np;
    msg << " header read" << endl;
    if(Np <= 0) {
        throw GeneralClassicException("GreenWakeFunction::setWakeFromFile",
                            " The particle number should be bigger than zero! Please check the first line of file \""
                            + filename_m +  "\".");
    }

    msg  << " Np = " << Np << endl;
    std::vector<double> wake(Np);
    std::vector<double> dist(Np);

    // read the wakefunction
    for(int i = 0; i < Np; i ++) {
        if(!fs.eof()) {
            fs >> dist[i] >> wake[i] >> dummy;
        }
        if(fs.eof()) {
            throw GeneralClassicException("GreenWakeFunction::setWakeFromFile",
                                " End of file reached before the whole wakefield is imported, please check file \""
                                + filename_m +  "\".");
        }
    }
    // if needed interpolate the wake in a way that the wake form the file fits to the wake needs in the code (??)

    FftWField_m.resize(NBin_m);

    for(int i = 0; i < NBin_m; i ++) {
        int j = 0;
        while(dist[j] < i * spacing) {
            j++;
        }
        -- j; // i * spacing should probably between dist[j] and dist[j+1]
        // linear interpolation
        FftWField_m[i] = wake[j] + ((wake[j+1] - wake[j]) / (dist[j+1] - dist[j]) * (i * spacing - dist[j]));
    }

    real = gsl_fft_real_wavetable_alloc(NBin_m);
    work = gsl_fft_real_workspace_alloc(NBin_m);

    gsl_fft_real_transform(FftWField_m.data(), 1, NBin_m, real, work);

    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
}

const std::string GreenWakeFunction::getType() const {
    return "GreenWakeFunction";
}
