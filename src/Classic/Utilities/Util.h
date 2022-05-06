//
// Namespace Util
//   This namespace contains useful global methods.
//
// Copyright (c) 200x - 2022, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#ifndef USEFULFUNCTIONS
#define USEFULFUNCTIONS

#include "Algorithms/Vektor.h"
#include "Algorithms/Quaternion.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <initializer_list>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>

// ------- DON'T DELETE: start --------
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define __DBGMSG__ __FILENAME__ << ": " << __LINE__ << "\t"
// ------- DON'T DELETE: end   --------

namespace Util {
    std::string getGitRevision();

    double erfinv(double x);

    inline
    double getGamma(Vector_t p) {
        return std::sqrt(dot(p, p) + 1.0);
    }

    inline
    Vector_t getBeta(Vector_t p) {
        return p / getGamma(p);
    }

    inline
    double getKineticEnergy(Vector_t p, double mass) {
        return (getGamma(p) - 1.0) * mass;
    }

    inline
    double getBetaGamma(double Ekin, double mass) {
        double value = std::sqrt(std::pow(Ekin / mass + 1.0, 2) - 1.0);
        if (value < std::numeric_limits<double>::epsilon())
            value = std::sqrt(2 * Ekin / mass);
        return value;
    }

    inline
    double convertMomentumEVoverCToBetaGamma(double p, double mass) {
        return p / mass;
    }

    inline
    std::string getTimeString(double time, unsigned int precision = 3) {
        std::string timeUnit(" [ps]");

        time *= 1e12;
        if (std::abs(time) > 1000) {
            time /= 1000;
            timeUnit = std::string(" [ns]");

            if (std::abs(time) > 1000) {
                time /= 1000;
                timeUnit = std::string(" [ms]");
            }
        } else if (std::abs(time) < 1.0) {
            time *= 1000;
            timeUnit = std::string(" [fs]");
        }

        std::stringstream timeOutput;
        timeOutput << std::fixed << std::setw(precision + 2) << std::setprecision(precision) << time << timeUnit;
        return timeOutput.str();
    }

    inline
    std::string getLengthString(double spos, unsigned int precision = 3) {
        std::string sposUnit(" [m]");

        if (std::abs(spos) < 1.0) {
            spos *= 1000.0;
            sposUnit = std::string(" [mm]");
        }

        if (std::abs(spos) < 1.0) {
            spos *= 1000.0;
            sposUnit = std::string(" [um]");
        }

        std::stringstream positionOutput;
        positionOutput << std::fixed << std::setw(precision + 2) << std::setprecision(precision) << spos << sposUnit;
        return positionOutput.str();
    }

    inline
    std::string getLengthString(Vector_t spos, unsigned int precision = 3) {
        std::string sposUnit(" [m]");
        double maxPos = std::abs(spos(0));
        for (unsigned int i = 1; i < 3u; ++ i) {
            maxPos = std::max(maxPos, std::abs(spos(i)));
        }

        std::stringstream positionOutput;

        if (maxPos < 1.0) {
            maxPos *= 1000.0;
            spos *= 1000.0;
            sposUnit = std::string(" [mm]");
        }

        if (maxPos < 1.0) {
            maxPos *= 1000.0;
            spos *= 1000.0;
            sposUnit = std::string(" [um]");
        }

        positionOutput << std::fixed << std::setprecision(precision)
                       << "( "
                       << std::setw(precision + 7) << spos(0) << " , "
                       << std::setw(precision + 7) << spos(1) << " , "
                       << std::setw(precision + 7) << spos(2)
                       << " )" << sposUnit;
        return positionOutput.str();
    }

    inline
    std::string getEnergyString(double energyInMeV, unsigned int precision = 3) {
        std::string energyUnit(" [MeV]");
        double energy = energyInMeV;

        if (energy > 1000.0) {
            energy /= 1000.0;
            energyUnit = std::string(" [GeV]");
        } else if (energy < 1.0) {
            energy *= 1000.0;
            energyUnit = std::string(" [keV]");
            if (energy < 1.0) {
                energy *= 1000.0;
                energyUnit = std::string(" [eV]");
            }
        }

        std::stringstream energyOutput;
        energyOutput << std::fixed << std::setw(precision + 2) << std::setprecision(precision) << energy << energyUnit;

        return energyOutput.str();
    }

    inline
    std::string getChargeString(double charge, unsigned int precision = 3) {
        std::string chargeUnit(" [fC]");

        charge *= 1e15;

        if (std::abs(charge) > 1000.0) {
            charge /= 1000.0;
            chargeUnit = std::string(" [pC]");
        }

        if (std::abs(charge) > 1000.0) {
            charge /= 1000.0;
            chargeUnit = std::string(" [nC]");
        }

        if (std::abs(charge) > 1000.0) {
            charge /= 1000.0;
            chargeUnit = std::string(" [uC]");
        }

        std::stringstream chargeOutput;
        chargeOutput << std::fixed << std::setw(precision + 2) << std::setprecision(precision) << charge << chargeUnit;

        return chargeOutput.str();
    }

    Vector_t getTaitBryantAngles(Quaternion rotation, const std::string& elementName = "");

    std::string toUpper(const std::string& str);

    std::string boolToUpperString(const bool& b);

    std::string boolVectorToUpperString(const std::vector<bool>& b);

    std::string doubleVectorToString(const std::vector<double>& v);

    std::string combineFilePath(std::initializer_list<std::string>);

    template<class IteratorIn, class IteratorOut>
    void toString(IteratorIn first, IteratorIn last, IteratorOut out);

    template <typename T>
    std::string toStringWithThousandSep(T value, char sep = '\'');

    struct KahanAccumulation {
        long double sum;
        long double correction;
        KahanAccumulation();

        KahanAccumulation& operator+=(double value);
    };

    unsigned int rewindLinesSDDS(const std::string& fileName, double maxSPos, bool checkForTime = true);

    std::string base64_encode(const std::string& string_to_encode);//unsigned char const* , unsigned int len);
    std::string base64_decode(std::string const& s);

    template<typename T, typename A>
    T* c_data(std::vector<T,A>& v) { return v.empty() ? static_cast<T*>(0) : &(v[0]); }

    template<typename T, typename A>
    T const* c_data(std::vector<T,A> const& v) { return v.empty() ? static_cast<T const*>(0) : &(v[0]); }
}

template <typename T>
std::string Util::toStringWithThousandSep(T value, char sep) {
    static_assert(std::is_integral<T>::value, "Util::toStringWithThousandSep: T must be of integer type");

    unsigned int powers = std::floor(std::max(0.0,
                                              std::log(std::abs((double)value)) / std::log(10.0))
                          );
    powers -= powers % 3u;

    std::ostringstream ret;
    unsigned int i = 0;
    while (powers >= 3u) {
        T multiplicator = std::pow(T(10), powers);
        T pre = value / multiplicator;
        if (i > 0) {
            ret << std::setw(3) << std::setfill('0') << pre << sep;
        } else {
            ret << pre << sep;
        }
        value -= pre * multiplicator;

        powers -= 3;
        ++ i;
    }

    if (i > 0) {
        ret << std::setw(3) << std::setfill('0') << value;
    } else {
        ret << value;
    }

    return ret.str();
}

template<class IteratorIn, class IteratorOut>
void Util::toString(IteratorIn first, IteratorIn last, IteratorOut out) {
    std::transform(first, last, out, [](auto d) {
        std::ostringstream stm;
        stm << d;
        return stm.str();
    } );
}

#endif
