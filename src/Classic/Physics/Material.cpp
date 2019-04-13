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

std::shared_ptr<Material> Material::addMaterial(const std::string &name,
                                                std::shared_ptr<Material> mat_ptr) {
    std::string nameUp = Util::toUpper(name);
    if (protoTable_sm.find(nameUp) != protoTable_sm.end())
        return protoTable_sm[nameUp];

    protoTable_sm.insert(std::make_pair(nameUp, mat_ptr));

    return mat_ptr;
}

std::shared_ptr<Material> Material::getMaterial(const std::string &name) {
    std::string nameUp = Util::toUpper(name);
    if (protoTable_sm.find(nameUp) != protoTable_sm.end()) return protoTable_sm[nameUp];

    throw GeneralClassicException("Material::getMaterial", "Unknown material '" + name + "'");
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
    auto berilium      = Material::addMaterial("Berilium",
                                               beryllium);
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
    auto titan         = Material::addMaterial("Titan",
                                               titanium);
    auto water         = Material::addMaterial("Water",
                                               std::shared_ptr<Material>(new Water()));
}