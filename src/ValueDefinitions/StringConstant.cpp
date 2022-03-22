//
// Class StringConstant
//   The STRING CONSTANT definition.
//
// Copyright (c) 2000 - 2021, Paul Scherrer Institut, Villigen PSI, Switzerland
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
#include "ValueDefinitions/StringConstant.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Utilities/Util.h"

#include <iostream>

#define CREATE_STRINGCONSTANT(x) \
    opal->create(new StringConstant(x, this, x));


StringConstant::StringConstant():
    ValueDefinition(1, "STRING_CONSTANT",
                    "The \"STRING CONSTANT\" statement defines a global "
                    "string constant:\n"
                    "\tSTRING CONSTANT <name> = <String-expression>;\n") {
    itsAttr[0] = Attributes::makeString("VALUE", "The constant value");

    registerOwnership(AttributeHandler::STATEMENT);

    OpalData *opal = OpalData::getInstance();
    opal->create(new StringConstant("GITREVISION", this, Util::getGitRevision()));

    // NOTE: Strings that contain a hyphen can't be written without quotes and
    // no string constant should be defined for them

    // Element / TYPE
    CREATE_STRINGCONSTANT("RING");
    CREATE_STRINGCONSTANT("CARBONCYCL");
    CREATE_STRINGCONSTANT("CYCIAE");
    CREATE_STRINGCONSTANT("AVFEQ");
    CREATE_STRINGCONSTANT("FFA");
    CREATE_STRINGCONSTANT("BANDRF");
    CREATE_STRINGCONSTANT("SYNCHROCYCLOTRON");
    CREATE_STRINGCONSTANT("SINGLEGAP");
    CREATE_STRINGCONSTANT("STANDING");
    CREATE_STRINGCONSTANT("TEMPORAL");
    CREATE_STRINGCONSTANT("SPATIAL");
    // Beam / PARTICLE
    CREATE_STRINGCONSTANT("ELECTRON");
    CREATE_STRINGCONSTANT("PROTON");
    CREATE_STRINGCONSTANT("POSITRON");
    CREATE_STRINGCONSTANT("ANTIPROTON");
    CREATE_STRINGCONSTANT("CARBON");
    CREATE_STRINGCONSTANT("HMINUS");
    CREATE_STRINGCONSTANT("URANIUM");
    CREATE_STRINGCONSTANT("MUON");
    CREATE_STRINGCONSTANT("DEUTERON");
    CREATE_STRINGCONSTANT("XENON");
    CREATE_STRINGCONSTANT("H2P");
    CREATE_STRINGCONSTANT("ALPHA");
    // Distribution / TYPE
    CREATE_STRINGCONSTANT("FROMFILE");
    CREATE_STRINGCONSTANT("GAUSS");
    CREATE_STRINGCONSTANT("BINOMIAL");
    CREATE_STRINGCONSTANT("FLATTOP");
    CREATE_STRINGCONSTANT("MULTIGAUSS");
    CREATE_STRINGCONSTANT("GUNGAUSSFLATTOPTH");
    CREATE_STRINGCONSTANT("ASTRAFLATTOPTH");
    CREATE_STRINGCONSTANT("GAUSSMATCHED");

    // Distribution / INPUTMOUNITS
    CREATE_STRINGCONSTANT("NONE");
    CREATE_STRINGCONSTANT("EVOVERC");

    // Distribution / EMISSIONMODEL
    CREATE_STRINGCONSTANT("ASTRA");
    CREATE_STRINGCONSTANT("NONEQUIL");
    // additionally: NONE

    // TrimCoil / TYPE
    // CREATE_STRINGCONSTANT("PSI-BFIELD");
    // CREATE_STRINGCONSTANT("PSI-PHASE");
    // CREATE_STRINGCONSTANT("PSI-BFIELD-MIRRORED");

    // OptimizeCmd / CROSSOVER
    CREATE_STRINGCONSTANT("BLEND");
    CREATE_STRINGCONSTANT("NAIVEONEPOINT");
    CREATE_STRINGCONSTANT("NAIVEUNIFORM");
    CREATE_STRINGCONSTANT("SIMULATEDBINARY");

    // OptimizeCmd / MUTATION
    CREATE_STRINGCONSTANT("ONEBIT");
    CREATE_STRINGCONSTANT("INDEPENDENTBIT");

    // OpalSample / TYPE
    CREATE_STRINGCONSTANT("UNIFORM");
    CREATE_STRINGCONSTANT("UNIFORM_INT");
    CREATE_STRINGCONSTANT("GAUSSIAN");
    CREATE_STRINGCONSTANT("LATIN_HYPERCUBE");

    // Option / RNGTYPE
    CREATE_STRINGCONSTANT("RANDOM");
    CREATE_STRINGCONSTANT("HALTON");
    CREATE_STRINGCONSTANT("SOBOL");
    CREATE_STRINGCONSTANT("NIEDERREITER");

    // OpalWake / TYPE
    // CREATE_STRINGCONSTANT("1D-CSR");
    // CREATE_STRINGCONSTANT("1D-CSR-IGF");
    // CREATE_STRINGCONSTANT("LONG-SHORT-RANGE");
    // CREATE_STRINGCONSTANT("TRANSV-SHORT-RANGE");

    // OpalWake / CONDUCT
    CREATE_STRINGCONSTANT("AC");
    CREATE_STRINGCONSTANT("DC");

    // ParticleMatterInteraction / TYPE
    CREATE_STRINGCONSTANT("SCATTERING");
    CREATE_STRINGCONSTANT("BEAMSTRIPPING");

    // FieldSolver / FSTYPE
    CREATE_STRINGCONSTANT("FFT");
    CREATE_STRINGCONSTANT("FFTPERIODIC");
    CREATE_STRINGCONSTANT("SAAMG");
    CREATE_STRINGCONSTANT("FMG");
    CREATE_STRINGCONSTANT("ML");
    CREATE_STRINGCONSTANT("AMR_MG");
    // additionally: NONE

    // FieldSolver / BCFFT
    CREATE_STRINGCONSTANT("OPEN");
    CREATE_STRINGCONSTANT("DIRICHLET");
    CREATE_STRINGCONSTANT("PERIODIC");

    // FieldSolver / GREENSF
    CREATE_STRINGCONSTANT("STANDARD");
    CREATE_STRINGCONSTANT("INTEGRATED");

    // FieldSolver / ITSOLVER
    CREATE_STRINGCONSTANT("CG");
    CREATE_STRINGCONSTANT("BICGSTAB");
    CREATE_STRINGCONSTANT("GMRES");

    // FieldSolver / INTERPL
    CREATE_STRINGCONSTANT("CONSTANT");
    CREATE_STRINGCONSTANT("LINEAR");
    CREATE_STRINGCONSTANT("QUADRATIC");

    // FieldSolver / PRECMODE
    CREATE_STRINGCONSTANT("STD");
    CREATE_STRINGCONSTANT("HIERARCHY");
    CREATE_STRINGCONSTANT("REUSE");

    // FieldSolver / AMR_MG_SMOOTHER
    CREATE_STRINGCONSTANT("GS");
    CREATE_STRINGCONSTANT("SGS");
    CREATE_STRINGCONSTANT("JACOBI");

    // FieldSolver / AMR_MG_PREC
    CREATE_STRINGCONSTANT("ILUT");
    CREATE_STRINGCONSTANT("CHEBYSHEV");
    CREATE_STRINGCONSTANT("RILUK");
    CREATE_STRINGCONSTANT("BLOCK_JACOBI");
    CREATE_STRINGCONSTANT("BLOCK_GS");
    CREATE_STRINGCONSTANT("SA");
    // additionally: NONE, JACOBI, GS

    // FieldSolver / AMR_MG_INTERP
    CREATE_STRINGCONSTANT("TRILINEAR");
    CREATE_STRINGCONSTANT("LAGRANGE");
    CREATE_STRINGCONSTANT("PC");

    // FieldSolver / AMR_MG_NORM
    CREATE_STRINGCONSTANT("L1_NORM");
    CREATE_STRINGCONSTANT("L2_NORM");
    CREATE_STRINGCONSTANT("LINF_NORM");

    // FieldSolver / AMR_MG_REUSE
    CREATE_STRINGCONSTANT("RP");
    CREATE_STRINGCONSTANT("RAP");
    // additionally: NONE

    // Option / PSDUMPFRAME
    CREATE_STRINGCONSTANT("GLOBAL");
    CREATE_STRINGCONSTANT("BUNCH_MEAN");
    CREATE_STRINGCONSTANT("REFERENCE");

    // OpalFilter / TYPE
    // CREATE_STRINGCONSTANT("SAVITZKY-GOLAY");
    CREATE_STRINGCONSTANT("FIXEDFFTLOWPASS");
    CREATE_STRINGCONSTANT("RELATIVEFFTLOWPASS");
    CREATE_STRINGCONSTANT("STENCIL");

    // TrackRun / METHOD
    CREATE_STRINGCONSTANT("THICK");
    // CREATE_STRINGCONSTANT("OPAL-T");
    // CREATE_STRINGCONSTANT("PARALLEL-T");
    // CREATE_STRINGCONSTANT("OPAL-CYCL");
    // CREATE_STRINGCONSTANT("CYCLOTRON-T");

    // TrackRun / MBMODE
    CREATE_STRINGCONSTANT("FORCE");
    CREATE_STRINGCONSTANT("AUTO");

    // TrackRun / MB_BINNING
    CREATE_STRINGCONSTANT("GAMMA_BINNING");
    CREATE_STRINGCONSTANT("BUNCH_BINNING");

    // ParticleMatterInteraction / MATERIAL
    CREATE_STRINGCONSTANT("AIR");
    CREATE_STRINGCONSTANT("ALUMINAAL2O3");
    CREATE_STRINGCONSTANT("ALUMINUM");
    CREATE_STRINGCONSTANT("BERYLLIUM");
    CREATE_STRINGCONSTANT("BORONCARBIDE");
    CREATE_STRINGCONSTANT("COPPER");
    CREATE_STRINGCONSTANT("GOLD");
    CREATE_STRINGCONSTANT("GRAPHITE");
    CREATE_STRINGCONSTANT("GRAPHITER6710");
    CREATE_STRINGCONSTANT("KAPTON");
    CREATE_STRINGCONSTANT("MOLYBDENUM");
    CREATE_STRINGCONSTANT("MYLAR");
    CREATE_STRINGCONSTANT("TITANIUM");
    CREATE_STRINGCONSTANT("WATER");

    // TrackCmd / TIMEINTEGRATOR
    CREATE_STRINGCONSTANT("RK_4");
    CREATE_STRINGCONSTANT("LF_2");
    CREATE_STRINGCONSTANT("MTS");

    // OpalVacuum / GAS
    CREATE_STRINGCONSTANT("H2");
    // additionally: AIR

    // BoundaryGeometry / TOPO
    CREATE_STRINGCONSTANT("RECTANGULAR");
    CREATE_STRINGCONSTANT("BOXCORNER");
    CREATE_STRINGCONSTANT("ELLIPTIC");

    // DumpEMFields / COORDINATE_SYSTEM
    CREATE_STRINGCONSTANT("CARTESIAN");
    CREATE_STRINGCONSTANT("CYLINDRICAL");
}

StringConstant::StringConstant(const std::string &name, StringConstant *parent):
    ValueDefinition(name, parent)
{}


StringConstant::StringConstant(const std::string &name, StringConstant *parent, const std::string &value):
    ValueDefinition(name, parent)
{
    Attributes::setString(itsAttr[0], value);
    itsAttr[0].setReadOnly(true);
    builtin = true;
}


StringConstant::~StringConstant()
{}


bool StringConstant::canReplaceBy(Object *) {
    return false;
}


StringConstant *StringConstant::clone(const std::string &name) {
    return new StringConstant(name, this);
}


std::string StringConstant::getString() const {
    return Attributes::getString(itsAttr[0]);
}


void StringConstant::print(std::ostream &os) const {
    os << "STRING " << getOpalName() << '=' << itsAttr[0] << ';';
    os << std::endl;
}

void StringConstant::printValue(std::ostream &os) const {
    os << itsAttr[0];
}