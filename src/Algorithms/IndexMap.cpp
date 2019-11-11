#include <map>
#include <limits>
#include <iostream>
#include <fstream>
#include <tuple>

#include "Algorithms/IndexMap.h"
#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Bend2D.h"
#include "Physics/Physics.h"
#include "Structure/ElementPositionWriter.h"

extern Inform *gmsg;

const double IndexMap::oneMinusEpsilon_m = 1.0 - std::numeric_limits<double>::epsilon();
namespace {
    void insertFlags(std::vector<double> &flags, std::shared_ptr<Component> element);
}

IndexMap::IndexMap():
    mapRange2Element_m(),
    mapElement2Range_m(),
    totalPathLength_m(0.0)
{ }

void IndexMap::print(std::ostream &out) const {
    if (mapRange2Element_m.size() == 0) return;

    out << "Size of map " << mapRange2Element_m.size() << " sections " << std::endl;
    out << std::fixed << std::setprecision(6);
    auto mapIti = mapRange2Element_m.begin();
    auto mapItf = mapRange2Element_m.end();

    double totalLength = (*mapRange2Element_m.rbegin()).first.second;
    unsigned int numDigits = std::floor(std::max(0.0, log(totalLength) / log(10.0))) + 1;

    for (; mapIti != mapItf; mapIti++) {
        const auto key = (*mapIti).first;
        const auto val = (*mapIti).second;
        out << "Key: ("
            << std::setw(numDigits + 7) << std::right << key.first
            << " - "
            << std::setw(numDigits + 7) << std::right << key.second
            << ") number of overlapping elements " << val.size() << "\n";

        for (auto element: val) {
            out << std::setw(25 + 2 * numDigits) << " " << element->getName() << "\n";
        }
    }
}

IndexMap::value_t IndexMap::query(key_t::first_type s, key_t::second_type ds) {
    const double lowerLimit = s - ds;//(ds < s? s - ds: 0);
    const double upperLimit = std::min(totalPathLength_m, s + ds);
    value_t elementSet;

    map_t::reverse_iterator rit = mapRange2Element_m.rbegin();
    if (rit != mapRange2Element_m.rend() && lowerLimit > (*rit).first.second) {
        throw OutOfBounds("IndexMap::query", "out of bounds");
    }

    map_t::iterator it = mapRange2Element_m.begin();
    const map_t::iterator end = mapRange2Element_m.end();

    for (; it != end; ++ it) {
        const double low = (*it).first.first;
        const double high = (*it).first.second;

        if (lowerLimit < high && upperLimit >= low) break;
    }

    if (it == end) return elementSet;

    map_t::iterator last = std::next(it);
    for (; last != end; ++ last) {
        const double low = (*last).first.first;

        if (upperLimit < low) break;
    }

    for (; it != last; ++ it) {
        const value_t &a = (*it).second;
        elementSet.insert(a.cbegin(), a.cend());
    }

    return elementSet;
}

void IndexMap::add(key_t::first_type initialS, key_t::second_type finalS, const value_t &val) {
    if (initialS > finalS) {
        std::swap(initialS, finalS);
    }
    key_t key(initialS, finalS * oneMinusEpsilon_m);

    mapRange2Element_m.insert(std::pair<key_t, value_t>(key, val));
    totalPathLength_m = (*mapRange2Element_m.rbegin()).first.second;

    value_t::iterator setIt = val.begin();
    const value_t::iterator setEnd = val.end();

    for (; setIt != setEnd; ++ setIt) {
        if (mapElement2Range_m.find(*setIt) == mapElement2Range_m.end()) {
            mapElement2Range_m.insert(std::make_pair(*setIt, key));
        } else {
            auto itpair = mapElement2Range_m.equal_range(*setIt);

            bool extendedExisting = false;
            for (auto it = itpair.first; it != itpair.second; ++ it) {
                key_t &currentRange = it->second;

                if (almostEqual(key.first, currentRange.second / oneMinusEpsilon_m)) {
                    currentRange.second = key.second;
                    extendedExisting = true;
                    break;
                }
            }
            if (!extendedExisting) {
                mapElement2Range_m.insert(std::make_pair(*setIt, key));
            }
        }
    }
}

void IndexMap::tidyUp(double zstop) {
    map_t::reverse_iterator rit = mapRange2Element_m.rbegin();

    if (rit != mapRange2Element_m.rend() &&
        (*rit).second.size() == 0 &&
        zstop > (*rit).first.first &&
        zstop < (*rit).first.second) {

        key_t key((*rit).first.first, zstop);
        value_t val;

        mapRange2Element_m.erase(std::next(rit).base());
        mapRange2Element_m.insert(std::pair<key_t, value_t>(key, val));
    }
}

enum elements {
    DIPOLE = 0,
    QUADRUPOLE,
    SEXTUPOLE,
    OCTUPOLE,
    DECAPOLE,
    MULTIPOLE,
    SOLENOID,
    RFCAVITY,
    MONITOR,
    OTHER,
    SIZE
};

void IndexMap::saveSDDS(double initialPathLength) const {
    if (mapRange2Element_m.size() == 0) return;

    std::vector<std::tuple<double, std::vector<double>, std::string> > sectors;

    // add for each sector four rows:
    // (s_i, 0)
    // (s_i, 1)
    // (s_f, 1)
    // (s_f, 0)
    // to the file, where
    // s_i is the start of the range and
    // s_f is the end of the range.
    auto mapIti = mapRange2Element_m.begin();
    auto mapItf = mapRange2Element_m.end();
    for (; mapIti != mapItf; mapIti++) {
        const auto &sectorElements = (*mapIti).second;
        if (sectorElements.size() == 0)
            continue;

        const auto &sectorRange = (*mapIti).first;

        double sectorBegin = sectorRange.first;
        double sectorEnd = sectorRange.second;

        std::vector<std::tuple<double, std::vector<double>, std::string> > currentSector(4);
        std::get<0>(currentSector[0]) = sectorBegin;
        std::get<0>(currentSector[1]) = sectorBegin;
        std::get<0>(currentSector[2]) = sectorEnd;
        std::get<0>(currentSector[3]) = sectorEnd;

        for (unsigned short i = 0; i < 4; ++ i) {
            auto &flags = std::get<1>(currentSector[i]);
            flags.resize(SIZE, 0);
        }

        for (auto element: sectorElements) {
            auto elementPassages = mapElement2Range_m.equal_range(element);
            auto passage = elementPassages.first;
            auto end = elementPassages.second;
            for (; passage != end; ++ passage) {
                const auto &elementRange = (*passage).second;
                double elementBegin = elementRange.first;
                double elementEnd = elementRange.second;

                if (elementBegin <= sectorBegin &&
                    elementEnd >= sectorEnd) {
                    break;
                }
            }

            const auto &elementRange = (*passage).second;
            if (elementRange.first < sectorBegin) {
                ::insertFlags(std::get<1>(currentSector[0]), element);
                std::get<2>(currentSector[0]) += element->getName() + ", ";
            }

            ::insertFlags(std::get<1>(currentSector[1]), element);
            std::get<2>(currentSector[1]) += element->getName() + ", ";

            ::insertFlags(std::get<1>(currentSector[2]), element);
            std::get<2>(currentSector[2]) += element->getName() + ", ";

            if (elementRange.second > sectorEnd) {
                ::insertFlags(std::get<1>(currentSector[3]), element);
                std::get<2>(currentSector[3]) += element->getName() + ", ";
            }
        }

        for (unsigned short i = 0; i < 4; ++ i) {
            sectors.push_back(currentSector[i]);
        }
    }

    // make the entries of the rf cavities a zigzag line
    const unsigned int numEntries = sectors.size();
    auto it = mapElement2Range_m.begin();
    auto end = mapElement2Range_m.end();
    for (; it != end; ++ it) {
        auto element = (*it).first;
        auto name = element->getName();
        auto type = element->getType();
        if (type != ElementBase::RFCAVITY &&
            type != ElementBase::TRAVELINGWAVE) {
            continue;
        }

        auto range = (*it).second;

        unsigned int i = 0;
        for (; i < numEntries; ++ i) {
            if (std::get<0>(sectors[i]) >= range.first) {
                break;
            }
        }

        if (i == numEntries) continue;

        unsigned int j = ++ i;
        while (std::get<0>(sectors[j]) < range.second) {
            ++ j;
        }

        double length = range.second - range.first;
        for (; i <= j; ++ i) {
            double pos = std::get<0>(sectors[i]);
            auto &items = std::get<1>(sectors[i]);

            items[RFCAVITY] = 1.0 - 2 * (pos - range.first) / length;
        }
    }

    // add row if range of first sector starts after initialPathLength
    if (sectors.size() > 0 &&
        std::get<0>(sectors[0]) > initialPathLength) {
        auto tmp = sectors;
        sectors = std::vector<std::tuple<double, std::vector<double>, std::string> >(1);
        std::get<0>(sectors[0]) = initialPathLength;
        std::get<1>(sectors[0]).resize(SIZE, 0.0);

        sectors.insert(sectors.end(), tmp.begin(), tmp.end());
    }

    std::string fileName("data/" + OpalData::getInstance()->getInputBasename() + "_ElementPositions.sdds");
    ElementPositionWriter writer(fileName);

    for (auto sector: sectors) {
        std::string names = std::get<2>(sector);
        if (names != "") {
            names = names.substr(0, names.length() - 2);
        }
        names = "\"" + names + "\"";
        writer.addRow(std::get<0>(sector),
                      std::get<1>(sector),
                      names);
    }
}

namespace {
    void insertFlags(std::vector<double> &flags, std::shared_ptr<Component> element) {
        switch (element->getType()) {
        case ElementBase::RBEND:
        case ElementBase::SBEND:
            {
                const Bend2D* bend = static_cast<const Bend2D*>(element.get());
                if (bend->getRotationAboutZ() > 0.5 * Physics::pi &&
                    bend->getRotationAboutZ() < 1.5 * Physics::pi) {
                    flags[DIPOLE] = -1;
                } else {
                    flags[DIPOLE] = 1;
                }
            }
            break;
        case ElementBase::MULTIPOLE:
            {
                const Multipole* mult = static_cast<const Multipole*>(element.get());
                switch(mult->getMaxNormalComponentIndex()) {
                case 1:
                    flags[DIPOLE] = (mult->isFocusing(0)? 1: -1);
                    break;
                case 2:
                    flags[QUADRUPOLE] = (mult->isFocusing(1)? 1: -1);
                    break;
                case 3:
                    flags[SEXTUPOLE] = (mult->isFocusing(2)? 1: -1);
                    break;
                case 4:
                    flags[OCTUPOLE] = (mult->isFocusing(3)? 1: -1);
                    break;
                case 5:
                    flags[DECAPOLE] = (mult->isFocusing(4)? 1: -1);
                    break;
                default:
                    flags[MULTIPOLE] = 1;
                }
            }
            break;
        case ElementBase::SOLENOID:
            flags[SOLENOID] = 1;
            break;
        case ElementBase::RFCAVITY:
        case ElementBase::TRAVELINGWAVE:
            flags[RFCAVITY] = 1;
            break;
        case ElementBase::MONITOR:
            flags[MONITOR] = 1;
            break;
        default:
            flags[OTHER] = 1;
            break;
        }

    }
}

std::pair<double, double> IndexMap::getRange(const IndexMap::value_t::value_type &element,
                                             double position) const {
    double minDistance = std::numeric_limits<double>::max();
    std::pair<double, double> range(0.0, 0.0);
    const std::pair<invertedMap_t::const_iterator, invertedMap_t::const_iterator> its = mapElement2Range_m.equal_range(element);
    if (std::distance(its.first, its.second) == 0)
        throw OpalException("IndexMap::getRange()",
                            "Element \"" + element->getName() + "\" not registered");

    for (invertedMap_t::const_iterator it = its.first; it != its.second; ++ it) {
        double distance = std::min(std::abs((*it).second.first - position),
                                   std::abs((*it).second.second - position));
        if (distance < minDistance) {
            minDistance = distance;
            range = (*it).second;
        }
    }

    return range;
}

IndexMap::value_t IndexMap::getTouchingElements(const std::pair<double, double> &range) {
    map_t::iterator it = mapRange2Element_m.begin();
    const map_t::iterator end = mapRange2Element_m.end();
    IndexMap::value_t touchingElements;

    for (; it != end; ++ it) {
        if (almostEqual(it->first.first, range.second) ||
            almostEqual(it->first.second, range.first))
            touchingElements.insert((it->second).begin(), (it->second).end());
    }

    return touchingElements;
}

bool IndexMap::almostEqual(double x, double y) {
    return (std::abs(x - y) < std::numeric_limits<double>::epsilon() * std::abs(x + y) * 2 ||
            std::abs(x - y) < std::numeric_limits<double>::min());
}