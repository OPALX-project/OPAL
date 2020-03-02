// ------------------------------------------------------------------------
// $RCSfile: OpalData.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.4.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalData
//   The global OPAL structure.
//   The OPAL object holds all global data required for a OPAL execution.
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:11 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "AbstractObjects/OpalData.h"
#include "Structure/DataSink.h"
#include "Structure/BoundaryGeometry.h"
#include "AbstractObjects/Attribute.h"
#include "AbstractObjects/Directory.h"
#include "AbstractObjects/Object.h"
#include "AbstractObjects/ObjectFunction.h"
#include "AbstractObjects/Table.h"
#include "AbstractObjects/ValueDefinition.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/RegularExpression.h"
#include <iostream>
#include <list>
#include <set>

// DTA
#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "ValueDefinitions/RealVariable.h"
#include "ValueDefinitions/StringConstant.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Parser/StringStream.h"
#include "Algorithms/PartBunchBase.h"
#include "Algorithms/bet/EnvelopeBunch.h"
// /DTA

#define MAX_NUM_INSTANCES 10
using namespace std;

// Class OpalData::ClearReference
// ------------------------------------------------------------------------

void OpalData::ClearReference::operator()(Object *object) const {
    object->clear();
}


// Struct OpalDataImpl.
// ------------------------------------------------------------------------

struct OpalDataImpl {
    OpalDataImpl();
    ~OpalDataImpl();

    // The main object directory.
    Directory mainDirectory;

    // The value of the global momentum.
    ValueDefinition *referenceMomentum;

    // The flag telling that something has changed.
    bool modified;

    // Directory of tables, for recalculation when something changes.
    std::list<Table *> tableDirectory;
    typedef std::list<Table *>::iterator tableIterator;

    // The set of expressions to be invalidated when something changes.
    std::set <AttributeBase *> exprDirectory;
    typedef std::set<AttributeBase *>::iterator exprIterator;

    // The page title from the latest TITLE command.
    std::string itsTitle_m;

    bool hasPriorRun_m;

    // true if we restart a simulation
    bool isRestart_m;

    // Input file name
    std::string inputFn_m;

    // Where to resume in a restart run
    int restartStep_m;

    // Where to resume in a restart run
    std::string restartFn_m;

    // true if the name of a restartFile is specified
    bool hasRestartFile_m;

    // dump frequency as found in restart file
    int restart_dump_freq_m;

    /// Mode for writing files
    OpalData::OPENMODE openMode_m = OpalData::OPENMODE::WRITE;

    // last step of a run
    int last_step_m;

    bool hasBunchAllocated_m;
    // The particle bunch to be tracked.
    PartBunchBase<double, 3> *bunch_m;

    bool hasDataSinkAllocated_m;

    DataSink *dataSink_m;

    bool hasSLBunchAllocated_m;
    // The particle bunch to be tracked.
    EnvelopeBunch *slbunch_m;

    // In units of seconds. This is half of tEmission
    double gPhaseShift_m;

    BoundaryGeometry *bg_m;

    std::vector<MaxPhasesT> maxPhases_m;
    energyEvolution_t energyEvolution_m;

    // The cartesian mesh
    Mesh_t *mesh_m;

    // The field layout f
    FieldLayout_t *FL_m;

    // The particle layout
    Layout_t *PL_m;

    // the accumulated (over all TRACKs) number of steps
    unsigned long long maxTrackSteps_m;

    bool isInOPALCyclMode_m;
    bool isInOPALTMode_m;
    bool isInOPALEnvMode_m;
    bool isOptimizerFlag_m;
    bool isInPrepState_m;

    std::map<std::string, unsigned int> problemSize_m;

    std::vector<std::string> arguments_m;
};


OpalDataImpl::OpalDataImpl():
    mainDirectory(), referenceMomentum(0), modified(false), itsTitle_m(),
    hasPriorRun_m(false),
    isRestart_m(false),
    restartStep_m(0),
    hasRestartFile_m(false),
    restart_dump_freq_m(1), last_step_m(0),
    hasBunchAllocated_m(false),
    hasDataSinkAllocated_m(false),
    hasSLBunchAllocated_m(false),
    gPhaseShift_m(0.0),
    maxTrackSteps_m(0),
    isInOPALCyclMode_m(false),
    isInOPALTMode_m(false),
    isInOPALEnvMode_m(false),
    isOptimizerFlag_m(false),
    isInPrepState_m(false)
{
    bunch_m    = nullptr;
    slbunch_m  = nullptr;
    dataSink_m = nullptr;
    bg_m       = nullptr;
    mesh_m     = nullptr;
    FL_m       = nullptr;
    PL_m       = nullptr;
}

OpalDataImpl::~OpalDataImpl() {
    // Make sure the main directory is cleared before the directories
    // for tables and expressions are deleted.
    delete mesh_m;// somehow this needs to be deleted first
    delete FL_m;
    //delete PL_m; //this gets deleted by FL_m

    delete bunch_m;
    delete slbunch_m;
    delete bg_m;
    delete dataSink_m;


    mainDirectory.erase();
    tableDirectory.clear();
    exprDirectory.clear();
}


// Class OpalData
// ------------------------------------------------------------------------

bool OpalData::isInstantiated = false;
OpalData *OpalData::instance = nullptr;
std::stack<OpalData*> OpalData::stashedInstances;

OpalData *OpalData::getInstance() {
    if(!isInstantiated) {
        instance = new OpalData();
        isInstantiated = true;
        return instance;
    } else {
        return instance;
    }
}

void OpalData::deleteInstance() {
    delete instance;
    instance = nullptr;
    isInstantiated = false;
}

void OpalData::stashInstance() {
    if (!isInstantiated) return;
    if (stashedInstances.size() + 1 > MAX_NUM_INSTANCES) {
        throw OpalException("OpalData::stashInstance()",
                            "too many OpalData instances stashed");
    }
    stashedInstances.push(instance);
    instance = nullptr;
    isInstantiated = false;
}

OpalData *OpalData::popInstance() {
    if (stashedInstances.size() == 0) {
        throw OpalException("OpalData::popInstance()",
                            "no OpalData instances stashed");
    }
    if (isInstantiated) deleteInstance();
    instance = stashedInstances.top();
    isInstantiated = true;
    stashedInstances.pop();

    return instance;
}

unsigned long long OpalData::getMaxTrackSteps() {
    return p->maxTrackSteps_m;
}

void OpalData::setMaxTrackSteps(unsigned long long s) {
    p->maxTrackSteps_m = s;
}

void OpalData::incMaxTrackSteps(unsigned long long s) {
    p->maxTrackSteps_m += s;
}


OpalData::OpalData() {
    p = new OpalDataImpl();
}

OpalData::~OpalData() {
    delete p;
}

void OpalData::reset() {
    delete p;
    p = new OpalDataImpl();
    p->hasPriorRun_m = false;
    p->isRestart_m = false;
    p->hasRestartFile_m = false;
    p->hasBunchAllocated_m = false;
    p->hasDataSinkAllocated_m = false;
    p->hasSLBunchAllocated_m = false;
    p->gPhaseShift_m = 0.0;
    p->maxPhases_m.clear();
    p->isInOPALCyclMode_m = false;
    p->isInOPALTMode_m = false;
    p->isInOPALEnvMode_m = false;
    p->isInPrepState_m = false;
    p->isOptimizerFlag_m = false;
}

bool OpalData::isInOPALCyclMode() {
    return p->isInOPALCyclMode_m;
}

bool OpalData::isInOPALTMode() {
    return  p->isInOPALTMode_m;
}

bool OpalData::isInOPALEnvMode() {
    return p->isInOPALEnvMode_m;
}

bool OpalData::isOptimizerRun() {
    return p->isOptimizerFlag_m;
}

void OpalData::setInOPALCyclMode() {
    p->isInOPALCyclMode_m = true;
}

void OpalData::setInOPALTMode() {
    p->isInOPALTMode_m = true;
}

void OpalData::setInOPALEnvMode() {
    p->isInOPALEnvMode_m = true;
}

void OpalData::setOptimizerFlag() {
    p->isOptimizerFlag_m = true;
}

bool OpalData::isInPrepState() {
    return p->isInPrepState_m;
}

void OpalData::setInPrepState(bool state) {
    p->isInPrepState_m = state;
}

bool OpalData::hasPriorTrack() {
    return p->hasPriorRun_m;
}

void OpalData::setPriorTrack(const bool &value) {
    p->hasPriorRun_m = value;
}

bool OpalData::inRestartRun() {
    return p->isRestart_m;
}

void OpalData::setRestartRun(const bool &value) {
    p->isRestart_m = value;
}


void OpalData::setRestartStep(int s) {
    p->restartStep_m = s;
}

int OpalData::getRestartStep() {
    return p->restartStep_m;
}


std::string OpalData::getRestartFileName() {
    return p->restartFn_m;
}


void OpalData::setRestartFileName(std::string s) {
    p->restartFn_m = s;
    p->hasRestartFile_m = true;
}

bool OpalData::hasRestartFile() {
    return p->hasRestartFile_m;

}

void OpalData::setRestartDumpFreq(const int &N) {
    p->restart_dump_freq_m = N;
}

int OpalData::getRestartDumpFreq() const {
    return p->restart_dump_freq_m;
}

void OpalData::setOpenMode(OPENMODE openMode) {
    p->openMode_m = openMode;
}

OpalData::OPENMODE OpalData::getOpenMode() const {
    return p->openMode_m;
}

void OpalData::setLastStep(const int &step) {
    p->last_step_m = step;
}

int OpalData::getLastStep() const {
    return p->last_step_m;
}

bool OpalData::hasSLBunchAllocated() {
    return p->hasSLBunchAllocated_m;
}

void OpalData::slbunchIsAllocated() {
    p->hasSLBunchAllocated_m = true;
}

void OpalData::setSLPartBunch(EnvelopeBunch *b) {
    p->slbunch_m = b;
}

EnvelopeBunch *OpalData::getSLPartBunch() {
    return p->slbunch_m;
}

bool OpalData::hasBunchAllocated() {
    return p->hasBunchAllocated_m;
}

void OpalData::bunchIsAllocated() {
    p->hasBunchAllocated_m = true;
}

void OpalData::setPartBunch(PartBunchBase<double, 3> *b) {
    p->bunch_m = b;
}

PartBunchBase<double, 3> *OpalData::getPartBunch() {
    return p->bunch_m;
}


bool OpalData::hasDataSinkAllocated() {
    return p->hasDataSinkAllocated_m;
}

void OpalData::setDataSink(DataSink *s) {
    p->dataSink_m = s;
    p->hasDataSinkAllocated_m = true;
}

DataSink *OpalData::getDataSink() {
    return p->dataSink_m;
}

void OpalData::setMaxPhase(std::string elName, double phi) {
    p->maxPhases_m.push_back(MaxPhasesT(elName, phi));
}

std::vector<MaxPhasesT>::iterator OpalData::getFirstMaxPhases() {
    return p->maxPhases_m.begin();
}

std::vector<MaxPhasesT>::iterator OpalData::getLastMaxPhases() {
    return p->maxPhases_m.end();
}

int OpalData::getNumberOfMaxPhases() {
    return p->maxPhases_m.size();
}

void OpalData::addEnergyData(double spos, double ekin) {
    p->energyEvolution_m.insert(std::make_pair(spos, ekin));
}

energyEvolution_t::iterator OpalData::getFirstEnergyData() {
    return p->energyEvolution_m.begin();
}

energyEvolution_t::iterator OpalData::getLastEnergyData() {
    return p->energyEvolution_m.end();
}


// Mesh_t* OpalData::getMesh() {
//  return p->mesh_m;
// }

// FieldLayout_t* OpalData::getFieldLayout() {
//  return p->FL_m;
// }

// Layout_t* OpalData::getLayout() {
//  return p->PL_m;
// }

// void OpalData::setMesh(Mesh_t *mesh) {
//  p->mesh_m = mesh;
// }

// void OpalData::setFieldLayout(FieldLayout_t *fieldlayout) {
//  p->FL_m = fieldlayout;
// }

// void OpalData::setLayout(Layout_t *layout) {
//  p->PL_m = layout;
// }

void OpalData::setGlobalPhaseShift(double shift) {
    /// units: (sec)
    p->gPhaseShift_m = shift;
}

double OpalData::getGlobalPhaseShift() {
    /// units: (sec)
    return p->gPhaseShift_m;
}

void   OpalData::setGlobalGeometry(BoundaryGeometry *bg) {
    p->bg_m = bg;
}

BoundaryGeometry *OpalData::getGlobalGeometry() {
    return p->bg_m;
}

bool OpalData::hasGlobalGeometry() {
    return p->bg_m != nullptr;
}



void OpalData::apply(const ObjectFunction &fun) {
    for(ObjectDir::iterator i = p->mainDirectory.begin();
        i != p->mainDirectory.end(); ++i) {
        fun(&*i->second);
    }
}


void OpalData::create(Object *newObject) {
    // Test for existing node with same name.
    const std::string name = newObject->getOpalName();
    Object *oldObject = p->mainDirectory.find(name);

    if(oldObject != nullptr) {
        throw OpalException("OpalData::create()",
                            "You cannot replace the object \"" + name + "\".");
    } else {
        p->mainDirectory.insert(name, newObject);
    }
}


void OpalData::define(Object *newObject) {
    // Test for existing node with same name.
    const std::string name = newObject->getOpalName();
    Object *oldObject = p->mainDirectory.find(name);

    if(oldObject != nullptr  &&  oldObject != newObject) {
        // Attempt to replace an object.
        if(oldObject->isBuiltin()  ||  ! oldObject->canReplaceBy(newObject)) {
            throw OpalException("OpalData::define()",
                                "You cannot replace the object \"" + name + "\".");
        } else {
            if(Options::info) {
                INFOMSG("Replacing the object \"" << name << "\"." << endl);
            }

            // Erase all tables which depend on the new object.
            OpalDataImpl::tableIterator i = p->tableDirectory.begin();
            while(i != p->tableDirectory.end()) {
                // We must increment i before calling erase(name),
                // since erase(name) removes "this" from "tables".
                Table *table = *i++;
                const std::string &tableName = table->getOpalName();

                if(table->isDependent(name)) {
                    if(Options::info) {
                        cerr << endl << "Erasing dependent table \"" << tableName
                             << "\"." << endl << endl;
                    }

                    // Remove table from directory.
                    // This erases the table from the main directory,
                    // and its destructor unregisters it from the table directory.
                    erase(tableName);
                }
            }

            // Replace all references to this object.
            for(ObjectDir::iterator i = p->mainDirectory.begin();
                i != p->mainDirectory.end(); ++i) {
                (*i).second->replace(oldObject, newObject);
            }

            // Remove old object.
            erase(name);
        }
    }

    // Force re-evaluation of expressions.
    p->modified = true;
    newObject->setDirty(true);
    p->mainDirectory.insert(name, newObject);

    // If this is a new definition of "P0", insert its definition.
    if(name == "P0") {
        if(ValueDefinition *p0 = dynamic_cast<ValueDefinition *>(newObject)) {
            setP0(p0);
        }
    }
}


void OpalData::erase(const std::string &name) {
    Object *oldObject = p->mainDirectory.find(name);

    if(oldObject != nullptr) {
        // Relink all children of "this" to "this->getParent()".
        for(ObjectDir::iterator i = p->mainDirectory.begin();
            i != p->mainDirectory.end(); ++i) {
            Object *child = &*i->second;
            if(child->getParent() == oldObject) {
                child->setParent(oldObject->getParent());
            }
        }

        // Remove the object.
        p->mainDirectory.erase(name);
    }
}


Object *OpalData::find(const std::string &name) {
    return p->mainDirectory.find(name);
}


double OpalData::getP0() const {
    static const double energy_scale = 1.0e+9;
    return p->referenceMomentum->getReal() * energy_scale;
}


void OpalData::makeDirty(Object *obj) {
    p->modified = true;
    if(obj) obj->setDirty(true);
}


void OpalData::printNames(std::ostream &os, const std::string &pattern) {
    int column = 0;
    RegularExpression regex(pattern);
    os << endl << "Object names matching the pattern \""
       << pattern << "\":" << endl;

    for(ObjectDir::const_iterator index = p->mainDirectory.begin();
        index != p->mainDirectory.end(); ++index) {
        const std::string name = (*index).first;

        if(! name.empty()  &&  regex.match(name)) {
            os << name;

            if(column < 80) {
                column += name.length();

                do {
                    os << ' ';
                    column++;
                } while((column % 20) != 0);
            } else {
                os << endl;
                column = 0;
            }
        }
    }

    if(column) os << endl;
    os << endl;
}


void OpalData::registerTable(Table *table) {
    p->tableDirectory.push_back(table);
}


void OpalData::unregisterTable(Table *table) {
    for(OpalDataImpl::tableIterator i = p->tableDirectory.begin();
        i != p->tableDirectory.end();) {
        OpalDataImpl::tableIterator j = i++;
        if(*j == table) p->tableDirectory.erase(j);
    }
}


void OpalData::registerExpression(AttributeBase *expr) {
    p->exprDirectory.insert(expr);
}


void OpalData::unregisterExpression(AttributeBase *expr) {
    p->exprDirectory.erase(expr);
}


void OpalData::setP0(ValueDefinition *p0) {
    p->referenceMomentum = p0;
}


void OpalData::storeTitle(const std::string &title) {
    p->itsTitle_m = title;
}

void OpalData::storeInputFn(const std::string &fn) {
    p->inputFn_m = fn;
}


void OpalData::printTitle(std::ostream &os) {
    os << p->itsTitle_m;
}

std::string OpalData::getTitle() {
    return p->itsTitle_m;
}

std::string OpalData::getInputFn() {
    return p->inputFn_m;
}

std::string OpalData::getInputBasename() {
    std::string & fn = p->inputFn_m;
    int const pdot = fn.rfind(".");
    return fn.substr(0, pdot);
}

void OpalData::update() {
    Inform msg("OpalData ");
    if(p->modified) {
        // Force re-evaluation of expressions.
        for(OpalDataImpl::exprIterator i = p->exprDirectory.begin();
            i != p->exprDirectory.end(); ++i) {
            (*i)->invalidate();
        }

        // Force refilling of dynamic tables.
        for(OpalDataImpl::tableIterator i = p->tableDirectory.begin();
            i != p->tableDirectory.end(); ++i) {
            (*i)->invalidate();
        }

        // Update all definitions.
        for(ObjectDir::iterator i = p->mainDirectory.begin();
            i != p->mainDirectory.end(); ++i) {
            (*i).second->update();
        }

        // Definitions are up-to-date.
        p->modified = false;
    }
}


std::map<std::string, std::string> OpalData::getVariableData() {
    std::map<std::string, std::string> udata;
    std::vector<std::string> uvars = this->getVariableNames();
    for (auto& uvar : uvars) {
        Object *tmpObject = OpalData::getInstance()->find(uvar);
        if (dynamic_cast<RealVariable*>(tmpObject)) {
            RealVariable* variable = dynamic_cast<RealVariable*>(OpalData::getInstance()->find(uvar));
            udata[uvar] = std::to_string(variable->getReal());
        } else if (dynamic_cast<StringConstant*>(tmpObject)) {
            StringConstant* variable = dynamic_cast<StringConstant*>(OpalData::getInstance()->find(uvar));
            udata[uvar] = variable->getString();
        } else {
            throw OpalException("OpalData::getVariableData()",
                                "Type of '" + uvar + "' not supported. "
                                "Only support for REAL and STRING.");
        }
    }
    return udata;
}

std::vector<std::string> OpalData::getVariableNames() {
    std::vector<std::string> result;

    for(ObjectDir::const_iterator index = p->mainDirectory.begin();
        index != p->mainDirectory.end(); ++index) {
        std::string tmpName = (*index).first;
        if (!tmpName.empty()) {
            Object *tmpObject = OpalData::getInstance()->find(tmpName);
            if (tmpObject) {
                if (!tmpObject || tmpObject->isBuiltin())
                    continue;
                if (tmpObject->getCategory() == "VARIABLE") {
                    result.push_back(tmpName);
                }
            }
        }
    }
    return result;
}

void OpalData::addProblemCharacteristicValue(const std::string &name, unsigned int value) {
    if (p->problemSize_m.find(name) != p->problemSize_m.end()) {
        p->problemSize_m.insert(std::make_pair(name, value));
    } else {
        p->problemSize_m[name] = value;
    }
}

const std::map<std::string, unsigned int> &OpalData::getProblemCharacteristicValues() const {
    return p->problemSize_m;
}

void OpalData::storeArguments(int argc, char *argv[]) {
    p->arguments_m.clear();
    for (int i = 0; i < argc; ++ i) {
        p->arguments_m.push_back(argv[i]);
    }
}

std::vector<std::string> OpalData::getArguments() {
    return p->arguments_m;
}