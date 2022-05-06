//
// Class Material
//   Base class for representing materials
//
// Copyright (c) 2019 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "Physics/Material.h"
#include "Physics/Air.h"
#include "Physics/AluminaAL2O3.h"
#include "Physics/Aluminum.h"
#include "Physics/Beryllium.h"
#include "Physics/BoronCarbide.h"
#include "Physics/Copper.h"
#include "Physics/Gold.h"
#include "Physics/Graphite.h"
#include "Physics/GraphiteR6710.h"
#include "Physics/Kapton.h"
#include "Physics/Molybdenum.h"
#include "Physics/Mylar.h"
#include "Physics/Titanium.h"
#include "Physics/Water.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"

#include <iostream>

using namespace Physics;

std::map<std::string, std::shared_ptr<Material> > Material::protoTable_sm;

std::shared_ptr<Material> Material::addMaterial(const std::string& name,
                                                std::shared_ptr<Material> mat_ptr) {
    std::string nameUp = Util::toUpper(name);
    if (protoTable_sm.find(nameUp) != protoTable_sm.end())
        return protoTable_sm[nameUp];

    protoTable_sm.insert(std::make_pair(nameUp, mat_ptr));

    return mat_ptr;
}

std::shared_ptr<Material> Material::getMaterial(const std::string& name) {
    std::string nameUp = Util::toUpper(name);
    return protoTable_sm[nameUp];
}

namespace {
    auto air           = Material::addMaterial("Air",
                                               std::shared_ptr<Material>(new Air()));
    auto aluminaal2o3  = Material::addMaterial("AluminaAL2O3",
                                               std::shared_ptr<Material>(new AluminaAL2O3()));
    auto aluminum      = Material::addMaterial("Aluminum",
                                               std::shared_ptr<Material>(new Aluminum()));
    auto beryllium     = Material::addMaterial("Beryllium",
                                               std::shared_ptr<Material>(new Beryllium()));
    auto boroncarbide  = Material::addMaterial("BoronCarbide",
                                               std::shared_ptr<Material>(new BoronCarbide()));
    auto copper        = Material::addMaterial("Copper",
                                               std::shared_ptr<Material>(new Copper()));
    auto gold          = Material::addMaterial("Gold",
                                               std::shared_ptr<Material>(new Gold()));
    auto graphite      = Material::addMaterial("Graphite",
                                               std::shared_ptr<Material>(new Graphite()));
    auto graphiter6710 = Material::addMaterial("GraphiteR6710",
                                               std::shared_ptr<Material>(new GraphiteR6710()));
    auto kapton        = Material::addMaterial("Kapton",
                                               std::shared_ptr<Material>(new Kapton()));
    auto molybdenum    = Material::addMaterial("Molybdenum",
                                               std::shared_ptr<Material>(new Molybdenum()));
    auto mylar         = Material::addMaterial("Mylar",
                                               std::shared_ptr<Material>(new Mylar()));
    auto titanium      = Material::addMaterial("Titanium",
                                               std::shared_ptr<Material>(new Titanium()));
    auto water         = Material::addMaterial("Water",
                                               std::shared_ptr<Material>(new Water()));
}