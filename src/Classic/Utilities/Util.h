#ifndef USEFULFUNCTIONS
#define USEFULFUNCTIONS

#include "Algorithms/Vektor.h"
#include "Algorithms/Quaternion.h"

#include <string>
#include <cstring>
#include <sstream>
#include <type_traits>
#include <functional>
#include <cmath>

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
    double getEnergy(Vector_t p, double mass) {
        return (getGamma(p) - 1.0) * mass;
    }

    inline
    double getP(double E, double mass) {
        double gamma = E / mass + 1;
        return std::sqrt(std::pow(gamma, 2.0) - 1.0);
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

    Vector_t getTaitBryantAngles(Quaternion rotation, const std::string &elementName = "");

    std::string toUpper(const std::string &str);

    template <typename T>
    std::string toStringWithThousandSep(T value, char sep = '\'');

    struct KahanAccumulation
    {
        long double sum;
        long double correction;
        KahanAccumulation();

        KahanAccumulation& operator+=(double value);
    };

    unsigned int rewindLinesSDDS(const std::string &fileName, double maxSPos, bool checkForTime = true);

    std::string base64_encode(const std::string &string_to_encode);//unsigned char const* , unsigned int len);
    std::string base64_decode(std::string const& s);
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

#endif