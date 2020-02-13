/*
 *  Copyright (c) 2014, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "Physics/Physics.h"
#include "Algorithms/PartBunchBase.h"
#include "AbsBeamline/BeamlineVisitor.h"

#include "AbsBeamline/EndFieldModel/EndFieldModel.h"
#include "AbsBeamline/VariableRFCavityFringeField.h"

VariableRFCavityFringeField::VariableRFCavityFringeField(const std::string &name) : VariableRFCavity(name) {
    initNull();  // initialise everything to NULL
}

VariableRFCavityFringeField::VariableRFCavityFringeField() : VariableRFCavity() {
    initNull();  // initialise everything to NULL
}

VariableRFCavityFringeField::VariableRFCavityFringeField(const VariableRFCavityFringeField& var) : VariableRFCavity() {
    initNull();  // initialise everything to NULL
    *this = var;
}

VariableRFCavityFringeField& VariableRFCavityFringeField::operator=(const VariableRFCavityFringeField& rhs) {
  
    if (&rhs == this) {
        return *this;
    }
    VariableRFCavity::operator=(rhs);
    setEndField(rhs.endField_m);
    zCentre_m = rhs.zCentre_m;
    f_m = rhs.f_m;
    g_m = rhs.f_m;
    h_m = rhs.f_m;
    return *this;
}

VariableRFCavityFringeField::~VariableRFCavityFringeField() {
}

ElementBase* VariableRFCavityFringeField::clone() const {
    return new VariableRFCavityFringeField(*this);
}

void VariableRFCavityFringeField::initNull() {
    VariableRFCavity::initNull();
    endField_m = std::shared_ptr<endfieldmodel::EndFieldModel>();
    zCentre_m = 0;
    maxOrder_m = 1;
}

void VariableRFCavityFringeField::accept(BeamlineVisitor& visitor) const {
    VariableRFCavity::initialise();
    VariableRFCavityFringeField* cavity =
                                 const_cast<VariableRFCavityFringeField*>(this);
    cavity->initialiseCoefficients();
    visitor.visitVariableRFCavity(*this);
}

bool VariableRFCavityFringeField::apply(const Vector_t &R, const Vector_t &/*P*/, const double &t, Vector_t &E, Vector_t &B) {
    if (R[2] > _length || R[2] < 0.) {
        return true;
    }
    if (std::abs(R[0]) > halfWidth_m || std::abs(R[1]) > halfHeight_m) {
        return true;
    }
    double z_pos = R[2]-zCentre_m;
    double E0 = amplitudeTD_m->getValue(t);
    double omega = Physics::two_pi*frequencyTD_m->getValue(t) * 1.0E-3; // need GHz on the element we have MHz
    double phi = phaseTD_m->getValue(t);
    double E_sin_t = E0*sin(omega * t + phi);
    double B_cos_t = E0*cos(omega * t + phi); // 1/c^2 factor in the h_n coefficients

    std::vector<double> y_power(maxOrder_m+1, 0.);
    y_power[0] = 1.;
    for (size_t i = 1; i < y_power.size(); ++i) {
        y_power[i] = y_power[i-1]*R[1];
    }

    // d^i f0 dz^i
    std::vector<double> endField(maxOrder_m/2+2, 0.);
    for (size_t i = 0; i < endField.size(); ++i) {
        endField[i] = endField_m->function(z_pos, i);
    }

    // omega^i
    std::vector<double> omegaPower(maxOrder_m+1, 1.);
    for (size_t i = 1; i < omegaPower.size(); ++i) {
        omegaPower[i] = omegaPower[i-1]*omega;
    }

    E = Vector_t(0., 0., 0.);
    B = Vector_t(0., 0., 0.);
    // even power of y
    //std::cerr << "EVEN POWER OF Y maxOrder: " << maxOrder_m << std::endl;
    for (size_t n = 0; n <= maxOrder_m ; n += 2) { // power of y
        double fCoeff = 0.;
        size_t index = n/2;
        //std::cerr << "Size i: " << index << " f_m[i]: " << f_m[index].size() 
        //          << " endfield: " << endField.size() << std::endl;
        for (size_t i = 0; i < f_m[index].size() && i < endField.size(); i += 2) { // derivative of f
            fCoeff += f_m[index][i]*endField[i]*omegaPower[n-i];
        }
        E[2] += E_sin_t*y_power[n]*fCoeff;
    }
    // odd power of y
    //std::cerr << "ODD POWER OF Y maxOrder: " << maxOrder_m << std::endl;
    for (size_t n = 1; n <= maxOrder_m; n += 2) {
        double gCoeff = 0.;
        double hCoeff = 0.;
        size_t index = (n-1)/2;
        //std::cerr << "Size i: " << index << " g_m[i]: " << g_m[index].size() << " endfield: " << endField.size() << std::endl;
        for (size_t j = 0; j < g_m[index].size() && j < endField.size(); ++j) {
            //std::cerr << "g_m        " << j << " " << g_m[index][j] << std::endl;
            //std::cerr << "endfield   " << j << " " << endField[j] << std::endl;
            //std::cerr << "omegaPower " << j << " " << omegaPower[n-j] << std::endl;
            gCoeff += g_m[index][j]*endField[j]*omegaPower[n-j];
        }
        for (size_t j = 0; j < h_m[index].size() && j < endField.size(); ++j) {
            hCoeff += h_m[index][j]*endField[j]*omegaPower[n-j];
            //std::cerr << "j: " << j << " " << hCoeff << " ";
        }
        //std::cerr << std::endl;
        E[1] += E_sin_t*y_power[n]*gCoeff;
        B[0] += B_cos_t*y_power[n]*hCoeff;
        //std::cerr << "APPLY B " << n << " " << B[0] << " " << hCoeff << std::endl;
    }
    B *= 1e1; //B converted to kGauss
    return false;
}

bool VariableRFCavityFringeField::apply(const size_t &i, const double &t,
                             Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], RefPartBunch_m->P[i], t, E, B);
}

bool VariableRFCavityFringeField::applyToReferenceParticle(const Vector_t &R,
                                                           const Vector_t &P,
                                                           const double &t,
                                                           Vector_t &E,
                                                           Vector_t &B) {
    return apply(R, P, t, E, B);
}

void VariableRFCavityFringeField::initialise(PartBunchBase<double, 3> *bunch,
                                             double &startField,
                                             double &endField) {
    VariableRFCavity::initialise(bunch, startField, endField);
    initialiseCoefficients();
}

void printVector(std::ostream& out, std::vector< std::vector<double> > vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        out << std::setw(3) << i;
        for (size_t j = 0; j < vec[i].size(); ++j) {
          out << " " << std::setw(14) << vec[i][j];
        }
        out << std::endl;
    }
}

void VariableRFCavityFringeField::initialiseCoefficients() {
    f_m = std::vector< std::vector<double> >();
    g_m = std::vector< std::vector<double> >();
    h_m = std::vector< std::vector<double> >();
    f_m.push_back(std::vector<double>(1, 1.));
    double c_l = Physics::c*1e-6; // speed of light in mm/ns
    // generate f_{n+2} term
    // note frequency term has to be added at apply(time) as we have
    // time-dependent frequency
    for(size_t n = 0; n+2 <= maxOrder_m; n += 2) {
        // n denotes the subscript on f_n
        // n+2 is the subscript on g_{n+2} and terms proportional to y^{n+2}
        std::vector<double> f_n = f_m.back(); // f_n
        std::vector<double> f_np2 = std::vector<double>(f_n.size()+2, 0.); // f_{n+2}
        double n_const = 1./(n+1.)/(n+2.);
        for (size_t j = 0; j < f_n.size(); ++j) {
            f_np2[j] -= f_n[j]*n_const/c_l/c_l;
        }
        for (size_t j = 0; j < f_n.size(); ++j) {
            f_np2[j+2] -= f_n[j]*n_const;
        }
        f_m.push_back(f_np2);
    }
    // generate g_{n+2} and h_{n+2} term
    for(size_t n = 0; n+1 <= maxOrder_m; n += 2) {
        // n denotes the subscript on f_n
        // n+1 is the subscript on g_{n+1} and terms proportional to y^{n+1}
        size_t f_index = n/2;
        std::vector<double> f_n = f_m[f_index];
        std::vector<double> g_np1 = std::vector<double>(f_n.size()+1, 0.);
        std::vector<double> h_np1 = std::vector<double>(f_n.size(), 0.);
        for (size_t j = 0; j < f_n.size(); ++j) {
            g_np1[j+1] = -1./(n+1.)*f_n[j];
            h_np1[j] = -1./c_l/c_l/(n+1.)*f_n[j];
        }
        g_m.push_back(g_np1);
        h_m.push_back(h_np1);
    }
}

void VariableRFCavityFringeField::printCoefficients(std::ostream& out) const {
    out << "f_m" << std::endl;
    printVector(out, f_m);
    out << "g_m" << std::endl;
    printVector(out, g_m);
    out << "h_m" << std::endl;
    printVector(out, h_m);
    out << std::endl;
}

void VariableRFCavityFringeField::setEndField(
                            std::shared_ptr<endfieldmodel::EndFieldModel> end) {
    endField_m = end;
}
