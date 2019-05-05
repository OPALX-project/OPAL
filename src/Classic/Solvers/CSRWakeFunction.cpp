#include "Solvers/CSRWakeFunction.hh"
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

CSRWakeFunction::CSRWakeFunction(const std::string &name, ElementBase *element, std::vector<Filter *> filters, const unsigned int &N):
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

void CSRWakeFunction::apply(PartBunchBase<double, 3> *bunch) {
    Inform msg("CSRWake ");

    const double sPos = bunch->get_sPos();
    std::pair<double, double> meshInfo;
    calculateLineDensity(bunch, meshInfo);
    const double &meshSpacing = meshInfo.second;
    const double &meshOrigin = meshInfo.first + 0.5 * meshSpacing;

    if(Ez_m.size() < lineDensity_m.size()) {
        Ez_m.resize(lineDensity_m.size(), 0.0);
        Psi_m.resize(lineDensity_m.size(), 0.0);
    }

    Vector_t smin, smax;
    bunch->get_bounds(smin, smax);
    double minPathLength = smin(2) + sPos - FieldBegin_m;
    if (sPos + smax(2) < FieldBegin_m) return;

    Ez_m[0] = 0.0;
    // calculate wake field of bunch
    for(unsigned int i = 1; i < lineDensity_m.size(); ++i) {
        Ez_m[i] = 0.0;

        double angleOfSlice = 0.0;
        double pathLengthOfSlice = minPathLength + i * meshSpacing;
        if(pathLengthOfSlice > 0.0)
            angleOfSlice = pathLengthOfSlice / bendRadius_m;

        calculateContributionInside(i, angleOfSlice, meshSpacing);
        calculateContributionAfter(i, angleOfSlice, meshSpacing);

        Ez_m[i] /= (4. * Physics::pi * Physics::epsilon_0);
    }

    // calculate the wake field seen by the particles
    for(unsigned int i = 0; i < bunch->getLocalNum(); ++i) {
        const Vector_t &R = bunch->R[i];
        double distanceToOrigin = (R(2) - meshOrigin) / meshSpacing;

        unsigned int indexz = (unsigned int)floor(distanceToOrigin);
        double leverz = distanceToOrigin - indexz;
        PAssert_LT(indexz, lineDensity_m.size() - 1);

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
            if (Ippl::myNode() == 0) {
                std::stringstream filename_str;
                filename_str << "data/" << bendName_m << "-CSRWake" << std::setw(5) << std::setfill('0') << file_number << ".txt";

                std::ofstream csr(filename_str.str().c_str());
                csr << std::setprecision(8);
                csr << "# " << sPos + smin(2) - FieldBegin_m << "\t" << sPos + smax(2) - FieldBegin_m << std::endl;
                for(unsigned int i = 0; i < lineDensity_m.size(); ++ i) {
                  csr << i *meshSpacing << "\t"
                      << Ez_m[i] << "\t"
                      << lineDensity_m[i] << "\t"
                      << dlineDensitydz_m[i] << std::endl;
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

void CSRWakeFunction::initialize(const ElementBase *ref) {
    double End;
    if(ref->getType() == ElementBase::RBEND ||
       ref->getType() == ElementBase::SBEND) {
        const Bend *bend = static_cast<const Bend *>(ref);
        bendRadius_m = bend->getBendRadius();
        bend->getDimensions(Begin_m, End);
        Length_m = bend->getEffectiveLength();
        FieldBegin_m = bend->getEffectiveCenter() - Length_m / 2.0;
        totalBendAngle_m = std::abs(bend->getBendAngle());
        bendName_m = bend->getName();
    }
}

void CSRWakeFunction::calculateLineDensity(PartBunchBase<double, 3> *bunch, std::pair<double, double> &meshInfo) {
    bunch->calcLineDensity(nBins_m, lineDensity_m, meshInfo);

    std::vector<Filter *>::const_iterator fit;
    for(fit = filters_m.begin(); fit != filters_m.end(); ++ fit) {
        (*fit)->apply(lineDensity_m);
    }

    dlineDensitydz_m.assign(lineDensity_m.begin(), lineDensity_m.end());
    diffOp_m->calc_derivative(dlineDensitydz_m, meshInfo.second);
}

void CSRWakeFunction::calculateContributionInside(size_t sliceNumber, double angleOfSlice, double meshSpacing) {
    if(angleOfSlice > totalBendAngle_m || angleOfSlice < 0.0) return;

    const double meshSpacingsup = pow(meshSpacing, -1. / 3.);
    double SlippageLength = pow(angleOfSlice, 3) * bendRadius_m / 24.;
    double relativeSlippageLength = SlippageLength / meshSpacing;
    if(relativeSlippageLength > sliceNumber) {

        /*
          Break integral into sum of integrals between grid points, then
          use linear interpolation between each grid point.
        */

        double dx1 = pow(sliceNumber, 2. / 3.);
        double dx2 = pow(sliceNumber, 5. / 3.);
        double dx3 = pow(sliceNumber - 1., 5. / 3.);
        Ez_m[sliceNumber] += 0.3 * meshSpacingsup * dlineDensitydz_m[0] * (5. * dx1 - 3. * dx2 + 3. * dx3);
        for(unsigned int j = 1; j < sliceNumber; ++ j) {
            dx1 = dx2;
            dx2 = dx3;
            dx3 = pow(sliceNumber - j - 1., 5. / 3.);
            Ez_m[sliceNumber] += 0.9 * meshSpacingsup * dlineDensitydz_m[j] * (dx1 - 2.* dx2 + dx3);
        }
        Ez_m[sliceNumber] += 0.9 * meshSpacingsup * dlineDensitydz_m[sliceNumber];

    } else if(relativeSlippageLength < 1) {

        // First do transient term.
        if(4.0 * relativeSlippageLength <= 1) {

            Ez_m[sliceNumber] += 3.0 * pow(SlippageLength, 2.0 / 3.0) * (lineDensity_m[sliceNumber] - lineDensity_m[sliceNumber - 1]) / meshSpacing;
        } else {

            if(4.0 * relativeSlippageLength < sliceNumber) {

                int j = sliceNumber - static_cast<int>(floor(4.0 * relativeSlippageLength));
                double frac = 4.0 * relativeSlippageLength - (sliceNumber - j);
                Ez_m[sliceNumber] -= (frac * lineDensity_m[j - 1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1. / 3.);

            }

            Ez_m[sliceNumber] += (relativeSlippageLength * lineDensity_m[sliceNumber - 1] + (1. - relativeSlippageLength) * lineDensity_m[sliceNumber]) / pow(SlippageLength, 1. / 3.);

        }

        // Now do steady state term.
        Ez_m[sliceNumber] += (0.3 / meshSpacing) * pow(SlippageLength, 2. / 3.) * (5. * dlineDensitydz_m[sliceNumber] - 2. * relativeSlippageLength * (dlineDensitydz_m[sliceNumber] - dlineDensitydz_m[sliceNumber - 1]));

    } else {

        if(4. * relativeSlippageLength < sliceNumber) {

            int j = sliceNumber - static_cast<int>(floor(4. * relativeSlippageLength));
            double frac = 4. * relativeSlippageLength - (sliceNumber - j);
            Ez_m[sliceNumber] -= (frac * lineDensity_m[j - 1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1. / 3.);

        }

        int j = sliceNumber - static_cast<int>(floor(SlippageLength / meshSpacing));
        double frac = relativeSlippageLength - (sliceNumber - j);
        Ez_m[sliceNumber] += (frac * lineDensity_m[j - 1] + (1. - frac) * lineDensity_m[j]) / pow(SlippageLength, 1. / 3.);

        double dx1 = pow(sliceNumber - j + frac, 2. / 3.);
        double dx2 = pow(sliceNumber - j, 2. / 3.);
        double dx3 = pow(sliceNumber - j + frac, 5. / 3.);
        double dx4 = pow(sliceNumber - j, 5. / 3.);

        Ez_m[sliceNumber] += 1.5 * meshSpacingsup * dlineDensitydz_m[j - 1] * (dx1 - dx2);
        Ez_m[sliceNumber] += 0.3 * meshSpacingsup * (dlineDensitydz_m[j] - dlineDensitydz_m[j - 1]) * (5.*(dx1 - dx2) + 3.*(dx3 - dx4));

        dx1 = dx2;
        dx2 = dx4;
        dx3 = pow(sliceNumber - j - 1., 5. / 3.);
        Ez_m[sliceNumber] += 0.3 * meshSpacingsup * dlineDensitydz_m[j] * (5.*dx1 - 3.*dx2 + 3.*dx3);
        for(unsigned int k = j + 1; k < sliceNumber; ++ k) {
            dx1 = dx2;
            dx2 = dx3;
            dx3 = pow(sliceNumber - k - 1., 5. / 3.);
            Ez_m[sliceNumber] += 0.9 * meshSpacingsup * dlineDensitydz_m[k] * (dx1 - 2.*dx2 + dx3);
        }
        Ez_m[sliceNumber] += 0.9 * meshSpacingsup * dlineDensitydz_m[sliceNumber];
    }
    double prefactor = -2. / pow(3. * bendRadius_m * bendRadius_m, 1. / 3.);
    Ez_m[sliceNumber] *= prefactor;
}

void CSRWakeFunction::calculateContributionAfter(size_t sliceNumber, double angleOfSlice, double meshSpacing) {
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

double CSRWakeFunction::calcPsi(const double &psiInitial, const double &x, const double &Ds) const {
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

const std::string CSRWakeFunction::getType() const {
    return "CSRWakeFunction";
}