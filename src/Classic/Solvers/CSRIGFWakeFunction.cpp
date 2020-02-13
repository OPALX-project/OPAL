#include "Solvers/CSRIGFWakeFunction.hh"
#include "Solvers/RootFinderForCSR.h"
#include "Algorithms/PartBunchBase.h"
#include "Filters/Filter.h"
#include "Filters/SavitzkyGolay.h"
#include "Physics/Physics.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/SBend.h"
#include "Utilities/Options.h"
#include "Utilities/Util.h"

#include <iostream>
#include <fstream>

extern Inform *gmsg;

CSRIGFWakeFunction::CSRIGFWakeFunction(const std::string &name, ElementBase *element, std::vector<Filter *> filters, const unsigned int &N):
    WakeFunction(name, element, N),
    filters_m(filters.begin(), filters.end()),
    lineDensity_m(),
    dlineDensitydz_m(),
    bendRadius_m(0.0),
    totalBendAngle_m(0.0)
{
    if (filters_m.size() == 0) {
        defaultFilter_m.reset(new SavitzkyGolayFilter(7, 3, 3, 3));
        filters_m.push_back(defaultFilter_m.get());
    }

    diffOp_m = filters_m.back();
}

void CSRIGFWakeFunction::apply(PartBunchBase<double, 3> *bunch) {
    Inform msg("CSRWake ");

    std::pair<double, double> meshInfo;
    calculateLineDensity(bunch, meshInfo);
    const double &meshOrigin = meshInfo.first;
    const double &meshSpacing = meshInfo.second;
    const unsigned int numOfSlices = lineDensity_m.size();

    if(Ez_m.size() < numOfSlices) {
        Ez_m.resize(numOfSlices, 0.0);
        Chi_m.resize(numOfSlices, 0.0);
        Grn_m.resize(numOfSlices, 0.0);
        Psi_m.resize(numOfSlices, 0.0);
    }

    for(unsigned int i = 0; i < numOfSlices; ++i) {
        Ez_m[i] = 0.0;
    }

    Vector_t smin, smax;
    bunch->get_bounds(smin, smax);
    double minPathLength = smin(2) + bunch->get_sPos() - FieldBegin_m;
    for(unsigned int i = 1; i < numOfSlices; i++) {
        double pathLengthOfSlice = minPathLength + i * meshSpacing;
        double angleOfSlice = pathLengthOfSlice/bendRadius_m;
        if (angleOfSlice > 0.0 && angleOfSlice <= totalBendAngle_m){
            calculateGreenFunction(bunch, meshSpacing);
        }
        // convolute with line density
        calculateContributionInside(i, angleOfSlice, meshSpacing);
        calculateContributionAfter(i, angleOfSlice, meshSpacing);
        Ez_m[i] /= 4 * Physics::pi * Physics::epsilon_0;
    }

    // calculate the wake field seen by the particles
    for(unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
        const Vector_t &R = bunch->R[i];
        unsigned int indexz = (unsigned int)floor((R(2) - meshOrigin) / meshSpacing);
        double leverz = (R(2) - meshOrigin) / meshSpacing - indexz;
        PAssert_LT(indexz + 1, numOfSlices);

        bunch->Ef[i](2) += (1. - leverz) * Ez_m[indexz] + leverz * Ez_m[indexz + 1];
    }

    if(Options::csrDump) {
        static std::string oldBendName;
        static unsigned long counter = 0;

        if(oldBendName != bendName_m) counter = 0;

        const int every = 1;
        bool print_criterion = (counter + 1) % every == 0;
        if(print_criterion) {
            static unsigned int file_number = 0;
            if(counter == 0) file_number = 0;
            double spos = bunch->get_sPos();
            if (Ippl::myNode() == 0) {
                std::stringstream filename_str;
                filename_str << "data/" << bendName_m << "-CSRWake" << std::setw(5) << std::setfill('0') << file_number << ".txt";
                std::ofstream csr(filename_str.str().c_str());
                csr << spos << ", " << FieldBegin_m << ", " << smin(2) << ", " << smax(2) << ", " << meshSpacing*64 << std::endl;
                for(unsigned int i = 0; i < lineDensity_m.size(); ++ i) {
                    csr << i *meshSpacing << "\t"
                        << Ez_m[i] << "\t"
                        << lineDensity_m[i] << std::endl;
                }
                csr.close();
                msg << "** wrote " << filename_str.str() << endl;
            }
            ++ file_number;
        }
        ++ counter;
        oldBendName = bendName_m;
    }
}

void CSRIGFWakeFunction::initialize(const ElementBase *ref) {
    if(ref->getType() == ElementBase::RBEND ||
       ref->getType() == ElementBase::SBEND) {

        const Bend2D *bend = static_cast<const Bend2D *>(ref);
        double End;

        bendRadius_m = bend->getBendRadius();
        bend->getDimensions(Begin_m, End);
        Length_m = bend->getEffectiveLength();
        FieldBegin_m = bend->getEffectiveCenter() - Length_m / 2.0;
        totalBendAngle_m = std::abs(bend->getBendAngle());
        bendName_m = bend->getName();
    }
}

void CSRIGFWakeFunction::calculateLineDensity(PartBunchBase<double, 3> *bunch, std::pair<double, double> &meshInfo) {
    bunch->calcLineDensity(nBins_m, lineDensity_m, meshInfo);

    // the following is only needed for after dipole
    std::vector<Filter *>::const_iterator fit;
    for(fit = filters_m.begin(); fit != filters_m.end(); ++ fit) {
        (*fit)->apply(lineDensity_m);
    }
    dlineDensitydz_m.assign(lineDensity_m.begin(), lineDensity_m.end());
    diffOp_m->calc_derivative(dlineDensitydz_m, meshInfo.second);
}

void CSRIGFWakeFunction::calculateGreenFunction(PartBunchBase<double, 3> *bunch, double meshSpacing)
{
    unsigned int numOfSlices = lineDensity_m.size();
    double gamma = bunch->get_meanKineticEnergy()/(bunch->getM()*1e-6)+1.0;
    double xmu_const = 3.0 * gamma * gamma * gamma / (2.0 * bendRadius_m);
    double chi_const = 9.0 / 16.0 * (6.0 - log(27.0 / 4.0));

    for(unsigned int i = 0; i < numOfSlices; ++i) {
        Chi_m[i] = 0.0;
        double z = i * meshSpacing;
        double xmu = xmu_const * z;
        double b = sqrt(xmu * xmu + 1.0) + xmu;
        if(xmu < 1e-3)
            Chi_m[i] = chi_const + 0.5 * pow(xmu, 2) - 7.0 / 54.0 * pow(xmu, 4) + 140.0 / 2187.0 * pow(xmu, 6);
        else
            Chi_m[i] = 9.0 / 16.0 * (3.0 * (-2.0 * xmu * pow(b, 1.0/3.0) + pow(b, 2.0/3.0) + pow(b, 4.0/3.0)) +
                                     log(pow((1 - pow(b, 2.0 / 3.0)) / xmu, 2) / (1 + pow(b, 2.0 / 3.0) + pow(b, 4.0 / 3.0))));
    }
    double grn_const = -16.0/(27.0 * gamma * gamma * meshSpacing);
    Grn_m[0] = grn_const * (Chi_m[1] - Chi_m[0]);
    Grn_m[numOfSlices - 1] = 0.0;
    for(unsigned int i = 1; i < numOfSlices - 1; ++i) {
        Grn_m[i] = grn_const * (Chi_m[i + 1] - 2.0 * Chi_m[i] + Chi_m[i - 1]);
    }
}

void CSRIGFWakeFunction::calculateContributionInside(size_t sliceNumber, double angleOfSlice, double /*meshSpacing*/)
{
    if(angleOfSlice > totalBendAngle_m || angleOfSlice < 0.0) return;
    int startSliceNum = 0;
    for(int j = sliceNumber; j >= startSliceNum; j--)
        Ez_m[sliceNumber] += lineDensity_m[j] * Grn_m[sliceNumber - j];
}

void CSRIGFWakeFunction::calculateContributionAfter(size_t sliceNumber, double angleOfSlice, double meshSpacing) {
    if(angleOfSlice <= totalBendAngle_m) return;

    double Ds_max = bendRadius_m * pow(totalBendAngle_m, 3) / 24. * (4. - 3.* totalBendAngle_m / angleOfSlice);

    // First do contribution from particles whose retarded position is
    // prior to the bend.
    double Ds_max2 = bendRadius_m * pow(totalBendAngle_m, 2) / 6. * (3. * angleOfSlice - 2. * totalBendAngle_m);
    int j = 0;
    double frac = 0.0;
    if(Ds_max2 / meshSpacing < sliceNumber) {
        j = sliceNumber - static_cast<int>(floor(Ds_max2 / meshSpacing));
        frac = Ds_max2 / meshSpacing - (sliceNumber - j);
        Ez_m[sliceNumber] -= (frac * lineDensity_m[j - 1] + (1. - frac) * lineDensity_m[j]) / (2. * angleOfSlice - totalBendAngle_m);
    }

    // Now do delta function contribution for particles whose retarded position
    // is in the bend.
    if(Ds_max / meshSpacing < sliceNumber) {
        j = sliceNumber - static_cast<int>(floor(Ds_max / meshSpacing));
        frac = Ds_max / meshSpacing - (sliceNumber - j);
        Ez_m[sliceNumber] += (frac * lineDensity_m[j - 1] + (1.0 - frac) * lineDensity_m[j]) / (2. * angleOfSlice - totalBendAngle_m);
    }

    // Now do integral contribution for particles whose retarded position is in
    // the bend.

    double angleOverlap = angleOfSlice - totalBendAngle_m;
    int k = sliceNumber;
    if(Ds_max / meshSpacing < sliceNumber) {
        k = j;
        Psi_m[k] = calcPsi(Psi_m[k], angleOverlap, meshSpacing * (k + frac));
        if(Psi_m[k] > 0 && Psi_m[k] < totalBendAngle_m)
            Ez_m[sliceNumber] += 0.5 * (frac * dlineDensitydz_m[sliceNumber - k - 1] + (1.0 - frac) * dlineDensitydz_m[sliceNumber - k]) / (Psi_m[k] + 2.0 * angleOverlap);
    } else {
        Psi_m[0] = calcPsi(Psi_m[0], angleOverlap, meshSpacing * sliceNumber);
        if(Psi_m[0] > 0 && Psi_m[0] < totalBendAngle_m)
            Ez_m[sliceNumber] += 0.5 * dlineDensitydz_m[0] / (Psi_m[0] + 2.0 * angleOverlap);
    }

    // Do rest of integral.
    for(unsigned int l = sliceNumber - k + 1; l < sliceNumber; ++ l) {
        Psi_m[l] = calcPsi(Psi_m[l], angleOverlap, meshSpacing * (sliceNumber - l));
        if(Psi_m[l] > 0 && Psi_m[l] < totalBendAngle_m)
            Ez_m[sliceNumber] += dlineDensitydz_m[l] / (Psi_m[l] + 2.0 * angleOverlap);
    }

    // We don't go right to the end as there is a singularity in the numerical integral that we don't quite know
    // how to deal with properly yet. This introduces a very slight error in the calculation (fractions of a percent).
    Psi_m[sliceNumber] = calcPsi(Psi_m[sliceNumber], angleOverlap, meshSpacing / 4.0);
    if(Psi_m[sliceNumber] > 0 && Psi_m[sliceNumber] < totalBendAngle_m)
        Ez_m[sliceNumber] += 0.5 * dlineDensitydz_m[sliceNumber] / (Psi_m[sliceNumber] + 2.0 * angleOverlap);

    double prefactor = -4 / bendRadius_m;
    Ez_m[sliceNumber] *= prefactor;
}

double CSRIGFWakeFunction::calcPsi(const double &psiInitial, const double &x, const double &Ds) const {
    /** solve the equation
     *  \f[
     *  \Delta s = \frac{R \Psi^3}{24} \frac{\Psi + 4x}{\Psi + x}
     *  \f]
     *  for \f$\Psi\f$ using Newtons method.
     */

    const int Nmax = 100;
    const double eps = 1e-10;
    double residual = 0.0;
    double psi = pow(24. * Ds / bendRadius_m, 1. / 3.);
    if(psiInitial != 0.0) psi = psiInitial;

    for(int i = 0; i < Nmax; ++i) {
        residual = bendRadius_m * psi * psi * psi * (psi + 4. * x) - 24. * Ds * psi - 24. * Ds * x;
        if(std::abs(residual) < eps)
            return psi;
        psi -= residual / (4. * bendRadius_m * psi * psi * psi + 12. * x * bendRadius_m * psi * psi - 24. * Ds);
    }
    RootFinderForCSR rootFinder(bendRadius_m, 4 * x * bendRadius_m, -24 * Ds, -24 * Ds * x);
    if (rootFinder.hasPositiveRealRoots()) {
        return rootFinder.searchRoot(eps);
    }

    ERRORMSG("In CSRWakeFunction::calcPsi(): exceed maximum number of iterations!" << endl);
    return psi;
}
const std::string CSRIGFWakeFunction::getType() const {
    return "CSRIGFWakeFunction";
}
