#include <map>
#include <limits>
#include <iostream>
#include <fstream>

#include "Algorithms/IndexMap.h"
#include "AbstractObjects/OpalData.h"
#include "AbsBeamline/Multipole.h"
#include "AbsBeamline/Bend2D.h"
#include "Physics/Physics.h"
#include "Structure/ElementPositionWriter.h"

extern Inform *gmsg;

const double IndexMap::oneMinusEpsilon_m = 1.0 - std::numeric_limits<double>::epsilon();

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
    const double lowerLimit = (ds < s? s - ds: 0);
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
    auto opal = OpalData::getInstance();
    if (opal->isOptimizerRun()) return;

    std::vector<std::vector<int> > allItems(SIZE);
    std::vector<double> allPositions;
    std::vector<std::string> allNames;
    unsigned int counter = 0;

    auto mapIti = mapRange2Element_m.begin();
    auto mapItf = mapRange2Element_m.end();
    for (; mapIti != mapItf; mapIti++) {
        const auto key = (*mapIti).first;
        const auto val = (*mapIti).second;

        std::vector<int> items(SIZE, 0);
        std::string names("");
        for (auto element: val) {
            switch (element->getType()) {
            case ElementBase::RBEND:
            case ElementBase::SBEND:
                {
                    const Bend2D* bend = static_cast<const Bend2D*>(element.get());
                    if (bend->getRotationAboutZ() > 0.5 * Physics::pi &&
                        bend->getRotationAboutZ() < 1.5 * Physics::pi) {
                        items[DIPOLE] = -1;
                    } else {
                        items[DIPOLE] = 1;
                    }
                }
                break;
            case ElementBase::MULTIPOLE:
                {
                    const Multipole* mult = static_cast<const Multipole*>(element.get());
                    switch(mult->getMaxNormalComponentIndex()) {
                    case 1:
                        items[DIPOLE] = (mult->isFocusing(0)? 1: -1);
                        break;
                    case 2:
                        items[QUADRUPOLE] = (mult->isFocusing(1)? 1: -1);
                        break;
                    case 3:
                        items[SEXTUPOLE] = (mult->isFocusing(2)? 1: -1);
                        break;
                    case 4:
                        items[OCTUPOLE] = (mult->isFocusing(3)? 1: -1);
                        break;
                    case 5:
                        items[DECAPOLE] = (mult->isFocusing(4)? 1: -1);
                        break;
                    default:
                        items[MULTIPOLE] = 1;
                    }
                }
                break;
            case ElementBase::SOLENOID:
                items[SOLENOID] = 1;
                break;
            case ElementBase::RFCAVITY:
            case ElementBase::TRAVELINGWAVE:
                items[RFCAVITY] = 1;
                break;
            case ElementBase::MONITOR:
                items[MONITOR] = 1;
                break;
            default:
                items[OTHER] = 1;
                break;
            }

            names += element->getName() + ", ";
        }

        if (names.length() == 0) continue;

        allPositions.push_back(key.first);
        allPositions.push_back(key.first);
        allPositions.push_back(key.second);
        allPositions.push_back(key.second);

        allNames.push_back("");
        allNames.push_back(names.substr(0, names.length() - 2));
        allNames.push_back(allNames.back());
        allNames.push_back("");

        for (unsigned int i = 0; i < SIZE; ++ i) {
            allItems[i].push_back(0);
            allItems[i].push_back(items[i]);
            allItems[i].push_back(items[i]);
            allItems[i].push_back(0);
        }

        ++ counter;
    }

    if (counter == 0) return;

    if (allPositions.front() > initialPathLength) {
        {
            auto tmp = allPositions;
            allPositions = std::vector<double>(4, initialPathLength);
            allPositions.insert(allPositions.end(), tmp.begin(), tmp.end());
        }

        {
            auto tmp = allNames;
            allNames = std::vector<std::string>(4, "");
            allNames.insert(allNames.end(), tmp.begin(), tmp.end());
        }

        for (unsigned int i = 0; i < SIZE; ++ i) {
            auto tmp = allItems[i];
            allItems[i] = std::vector<int>(4, 0);
            allItems[i].insert(allItems[i].end(), tmp.begin(), tmp.end());
        }
    }

    const unsigned int totalSize = counter;
    for (unsigned int i = 0; i < SIZE; ++ i) {
        for (unsigned int j = 0; j < totalSize - 1; ++ j) {
            if (allItems[i][4 * j + 1] == 0) continue;
            if (allItems[i][4 * (j + 1) + 1] == 0) continue;
            if (allItems[i][4 * j + 1] * allItems[i][4 * (j + 1) + 1] < 0) continue;
            if (std::abs(allPositions[4 * j + 3] - allPositions[4 * (j + 1)]) > 1e-6) continue;

            allItems[i][4 * j + 3] = allItems[i][4 * (j + 1)] = allItems[i][4 * j + 1];
        }
    }

    for (unsigned int j = 0; j < totalSize; ++ j) {
        allItems[RFCAVITY][4 * j + 2] = -allItems[RFCAVITY][4 * j + 1];
    }

    std::string fileName("data/" + OpalData::getInstance()->getInputBasename() + "_ElementPositions.sdds");
    ElementPositionWriter writer(fileName);

    for (unsigned int i = 0; i < 4 * totalSize; ++ i) {
        std::vector<int> values(SIZE, 0);
        for (unsigned int j = 0; j < SIZE; ++ j) {
            values[j] = allItems[j][i];
        }
        writer.addRow(allPositions[i],
                      values,
                      allNames[i]);
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