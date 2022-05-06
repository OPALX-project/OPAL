//
// Class IndexMap
//
// This class stores and prints the sequence of elements that the referenc particle passes.
// Each time the reference particle enters or leaves an element an entry is added to the map.
// With help of this map one can determine which element can be found at a given position.
//
// Copyright (c) 2016,       Christof Metzger-Kraus, Helmholtz-Zentrum Berlin, Germany
//               2017 - 2020 Christof Metzger-Kraus
//
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
#include "Utilities/Util.h"

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
    if (mapRange2Element_m.empty()) return;

    out << "* Size of map " << mapRange2Element_m.size() << " sections " << std::endl;
    out << std::fixed << std::setprecision(6);
    auto mapIti = mapRange2Element_m.begin();
    auto mapItf = mapRange2Element_m.end();

    double totalLength = (*mapRange2Element_m.rbegin()).first.end;
    unsigned int numDigits = std::floor(std::max(0.0, log(totalLength) / log(10.0))) + 1;

    for (; mapIti != mapItf; mapIti++) {
        const key_t key = (*mapIti).first;
        const value_t val = (*mapIti).second;
        out << "* Key: ("
            << std::setw(numDigits + 7) << std::right << key.begin
            << " - "
            << std::setw(numDigits + 7) << std::right << key.end
            << ") number of overlapping elements " << val.size() << "\n";

        for (auto element: val) {
            out << "* " << std::setw(25 + 2 * numDigits) << " " << element->getName() << "\n";
        }
    }
}

IndexMap::value_t IndexMap::query(key_t::first_type s, key_t::second_type ds) {
    const double lowerLimit = s - ds;//(ds < s? s - ds: 0);
    const double upperLimit = std::min(totalPathLength_m, s + ds);
    value_t elementSet;

    map_t::reverse_iterator rit = mapRange2Element_m.rbegin();
    if (rit != mapRange2Element_m.rend() && lowerLimit > (*rit).first.end) {
        throw OutOfBounds("IndexMap::query", "out of bounds");
    }

    map_t::iterator it = mapRange2Element_m.begin();
    const map_t::iterator end = mapRange2Element_m.end();

    for (; it != end; ++ it) {
        const double low = (*it).first.begin;
        const double high = (*it).first.end;

        if (lowerLimit < high && upperLimit >= low) break;
    }

    if (it == end) return elementSet;

    map_t::iterator last = std::next(it);
    for (; last != end; ++ last) {
        const double low = (*last).first.begin;

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
    key_t key{initialS, finalS * oneMinusEpsilon_m};

    mapRange2Element_m.insert(std::pair<key_t, value_t>(key, val));
    totalPathLength_m = (*mapRange2Element_m.rbegin()).first.end;

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

                if (almostEqual(key.begin, currentRange.end / oneMinusEpsilon_m)) {
                    currentRange.end = key.end;
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
        (*rit).second.empty() &&
        zstop > (*rit).first.begin) {

        key_t key{(*rit).first.begin, zstop};
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
    if (mapRange2Element_m.empty()) return;

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
        if (sectorElements.empty())
            continue;

        const auto &sectorRange = (*mapIti).first;

        double sectorBegin = sectorRange.begin;
        double sectorEnd = sectorRange.end;

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
                double elementBegin = elementRange.begin;
                double elementEnd = elementRange.end;

                if (elementBegin <= sectorBegin &&
                    elementEnd >= sectorEnd) {
                    break;
                }
            }

            const auto &elementRange = (*passage).second;
            if (elementRange.begin < sectorBegin) {
                ::insertFlags(std::get<1>(currentSector[0]), element);
                std::get<2>(currentSector[0]) += element->getName() + ", ";
            }

            ::insertFlags(std::get<1>(currentSector[1]), element);
            std::get<2>(currentSector[1]) += element->getName() + ", ";

            ::insertFlags(std::get<1>(currentSector[2]), element);
            std::get<2>(currentSector[2]) += element->getName() + ", ";

            if (elementRange.end > sectorEnd) {
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
        if (type != ElementType::RFCAVITY &&
            type != ElementType::TRAVELINGWAVE) {
            continue;
        }

        auto range = (*it).second;

        unsigned int i = 0;
        for (; i < numEntries; ++ i) {
            if (std::get<0>(sectors[i]) >= range.begin) {
                break;
            }
        }

        if (i == numEntries) continue;

        unsigned int j = ++ i;
        while (std::get<0>(sectors[j]) < range.end) {
            ++ j;
        }

        double length = range.end - range.begin;
        for (; i <= j; ++ i) {
            double pos = std::get<0>(sectors[i]);
            auto &items = std::get<1>(sectors[i]);

            items[RFCAVITY] = 1.0 - 2 * (pos - range.begin) / length;
        }
    }

    // add row if range of first sector starts after initialPathLength
    if (!sectors.empty() &&
        std::get<0>(sectors[0]) > initialPathLength) {
        auto tmp = sectors;
        sectors = std::vector<std::tuple<double, std::vector<double>, std::string> >(1);
        std::get<0>(sectors[0]) = initialPathLength;
        std::get<1>(sectors[0]).resize(SIZE, 0.0);

        sectors.insert(sectors.end(), tmp.begin(), tmp.end());
    }

    std::string fileName = Util::combineFilePath({
        OpalData::getInstance()->getAuxiliaryOutputDirectory(),
        OpalData::getInstance()->getInputBasename() + "_ElementPositions.sdds"
    });
    ElementPositionWriter writer(fileName);

    for (auto sector: sectors) {
        std::string names = std::get<2>(sector);
        if (!names.empty()) {
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
        case ElementType::RBEND:
        case ElementType::SBEND:
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
        case ElementType::MULTIPOLE:
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
        case ElementType::SOLENOID:
            flags[SOLENOID] = 1;
            break;
        case ElementType::RFCAVITY:
        case ElementType::TRAVELINGWAVE:
            flags[RFCAVITY] = 1;
            break;
        case ElementType::MONITOR:
            flags[MONITOR] = 1;
            break;
        default:
            flags[OTHER] = 1;
            break;
        }

    }
}

IndexMap::key_t IndexMap::getRange(const IndexMap::value_t::value_type &element,
                                   double position) const {
    double minDistance = std::numeric_limits<double>::max();
    key_t range{0.0, 0.0};
    const std::pair<invertedMap_t::const_iterator, invertedMap_t::const_iterator> its = mapElement2Range_m.equal_range(element);
    if (std::distance(its.first, its.second) == 0)
        throw OpalException("IndexMap::getRange()",
                            "Element \"" + element->getName() + "\" not registered");

    for (invertedMap_t::const_iterator it = its.first; it != its.second; ++ it) {
        double distance = std::min(std::abs((*it).second.begin - position),
                                   std::abs((*it).second.end - position));
        if (distance < minDistance) {
            minDistance = distance;
            range = (*it).second;
        }
    }

    return range;
}

IndexMap::value_t IndexMap::getTouchingElements(const IndexMap::key_t &range) const {
    map_t::const_iterator it = mapRange2Element_m.begin();
    const map_t::const_iterator end = mapRange2Element_m.end();
    value_t touchingElements;

    for (; it != end; ++ it) {
        if (almostEqual(it->first.begin, range.begin) ||
            almostEqual(it->first.end, range.end))
            touchingElements.insert((it->second).begin(), (it->second).end());
    }

    return touchingElements;
}

bool IndexMap::almostEqual(double x, double y) {
    return (std::abs(x - y) < std::numeric_limits<double>::epsilon() * std::abs(x + y) * 2 ||
            std::abs(x - y) < std::numeric_limits<double>::min());
}
