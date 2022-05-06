//
// Class OpalSample
//   The SAMPLING definition.
//   A SAMPLING definition is used to run the optimizer in sample mode.
//
// Copyright (c) 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
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
#include "Sample/OpalSample.h"

#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Sample/FromFile.h"
#include "Sample/LatinHyperCube.h"
#include "Sample/Normal.h"
#include "Sample/SampleGaussianSequence.h"
#include "Sample/SampleRandomizedSequence.h"
#include "Sample/SampleSequence.h"
#include "Sample/Uniform.h"
#include "Utilities/OpalException.h"
#include "Utilities/Util.h"

#include <unordered_map>

// The attributes of class OpalSample.
namespace {
    enum {
        TYPE,       // The type of sampling
        VARIABLE,   // name of design variable
        SEED,       // for random sample methods
        FNAME,      // file to read from sampling points
        N,
        RANDOM,
        STEP,
        SIZE
    };
}

OpalSample::OpalSample():
    Definition(SIZE, "SAMPLING",
               "The \"SAMPLING\" statement defines methods used for the optimizer in sample mode.")
    , size_m(1)
{
    itsAttr[TYPE] = Attributes::makePredefinedString
        ("TYPE", "Distribution type.",
         {"UNIFORM_INT",
          "UNIFORM",
          "GAUSSIAN",
          "FROMFILE",
          "LATIN_HYPERCUBE",
          "RANDOM_SEQUENCE_UNIFORM_INT",
          "RANDOM_SEQUENCE_UNIFORM"});

    itsAttr[VARIABLE] = Attributes::makeString
        ("VARIABLE", "Name of design variable");

    itsAttr[SEED] = Attributes::makeReal
        ("SEED", "seed for random sampling");

    itsAttr[FNAME] = Attributes::makeString
        ("FNAME", "File to read from the sampling points");

    itsAttr[N] = Attributes::makeReal
        ("N", "Number of sampling points", 1);

    itsAttr[RANDOM] = Attributes::makeBool
        ("RANDOM", "Whether sequence should be sampled randomly (default: false)", false);

    itsAttr[STEP] = Attributes::makeReal
        ("STEP", "Increment for randomized sequences (default: 1)", 1.0);

    registerOwnership(AttributeHandler::STATEMENT);
}


OpalSample::OpalSample(const std::string& name, OpalSample* parent):
    Definition(name, parent)
{}


OpalSample* OpalSample::clone(const std::string& name) {
    return new OpalSample(name, this);
}


void OpalSample::execute() {

}


OpalSample* OpalSample::find(const std::string& name) {
    OpalSample* sampling = dynamic_cast<OpalSample*>(OpalData::getInstance()->find(name));

    if (sampling == nullptr) {
        throw OpalException("OpalSample::find()",
                            "OpalSample \"" + name + "\" not found.");
    }
    return sampling;
}


void OpalSample::initialize(const std::string& dvarName,
                            double lower, double upper,
                            size_t modulo, bool /*sequence*/) {

    if ( lower >= upper ) {
        throw OpalException("OpalSample::initialize()",
                            "Lower bound >= upper bound.");
    }

    static const std::unordered_map<std::string, OpalSampleMethod> stringOpalSampleMethod_s = {
        {"UNIFORM_INT",                 OpalSampleMethod::UNIFORM_INT},
        {"UNIFORM",                     OpalSampleMethod::UNIFORM},
        {"GAUSSIAN",                    OpalSampleMethod::GAUSSIAN},
        {"FROMFILE",                    OpalSampleMethod::FROMFILE},
        {"LATIN_HYPERCUBE",             OpalSampleMethod::LATIN_HYPERCUBE},
        {"RANDOM_SEQUENCE_UNIFORM_INT", OpalSampleMethod::RANDOM_SEQUENCE_UNIFORM_INT},
        {"RANDOM_SEQUENCE_UNIFORM",     OpalSampleMethod::RANDOM_SEQUENCE_UNIFORM}
    };
    std::string type = Attributes::getString(itsAttr[TYPE]);
    if (type.empty()) {
        throw OpalException("OpalSample::initialize",
                            "The attribute \"TYPE\" isn't set for the \"SAMPLING\" statement");
    }
    OpalSampleMethod method = stringOpalSampleMethod_s.at(type);

    int seed = Attributes::getReal(itsAttr[SEED]);
    size_m = Attributes::getReal(itsAttr[N]);
    double step = Attributes::getReal(itsAttr[STEP]);
    bool random = Attributes::getBool(itsAttr[RANDOM]);

    if (!random) {
        if (method == OpalSampleMethod::UNIFORM_INT) {
            sampleMethod_m.reset( new SampleSequence<int>(lower, upper, modulo, size_m) );
        } else if (method == OpalSampleMethod::UNIFORM) {
            sampleMethod_m.reset( new SampleSequence<double>(lower, upper, modulo, size_m) );
        } else if (method == OpalSampleMethod::GAUSSIAN) {
            sampleMethod_m.reset( new SampleGaussianSequence(lower, upper, modulo, size_m) );
        } else if (method == OpalSampleMethod::FROMFILE) {
            std::string fname = Attributes::getString(itsAttr[FNAME]);
            sampleMethod_m.reset( new FromFile(fname, dvarName, modulo) );
        } else {
            throw OpalException("OpalSample::initialize",
                                "The sampling method \"TYPE=" + type + "\" is not supported out of random sampling mode");
        }
    } else {
        switch (method) {
            case OpalSampleMethod::UNIFORM_INT: {
                if (Attributes::getReal(itsAttr[SEED])) {
                    sampleMethod_m.reset( new Uniform<int>(lower, upper, seed) );
                } else {
                    sampleMethod_m.reset( new Uniform<int>(lower, upper) );
                }
                break;
            }
            case OpalSampleMethod::UNIFORM: {
                if (Attributes::getReal(itsAttr[SEED])) {
                    sampleMethod_m.reset( new Uniform<double>(lower, upper, seed) );
                } else {
                    sampleMethod_m.reset( new Uniform<double>(lower, upper) );
                }
                break;
            }
            case OpalSampleMethod::GAUSSIAN: {
                if (Attributes::getReal(itsAttr[SEED])) {
                    sampleMethod_m.reset( new Normal(lower, upper, seed) );
                } else {
                    sampleMethod_m.reset( new Normal(lower, upper) );
                }
                break;
            }
            case OpalSampleMethod::FROMFILE: {
                std::string fname = Attributes::getString(itsAttr[FNAME]);
                sampleMethod_m.reset( new FromFile(fname, dvarName, modulo) );
                size_m = static_cast<FromFile*>(sampleMethod_m.get())->getSize();
                break;
            }
            case OpalSampleMethod::LATIN_HYPERCUBE: {
                if (Attributes::getReal(itsAttr[SEED])) {
                    sampleMethod_m.reset( new LatinHyperCube(lower, upper, seed) );
                } else {
                    sampleMethod_m.reset( new LatinHyperCube(lower, upper) );
                }
                break;
            }
            case OpalSampleMethod::RANDOM_SEQUENCE_UNIFORM_INT: {
                if (Attributes::getReal(itsAttr[SEED])) {
                    sampleMethod_m.reset(
                        new SampleRandomizedSequence<int>(lower, upper, step, seed)
                    );
                } else {
                    sampleMethod_m.reset(
                        new SampleRandomizedSequence<int>(lower, upper, step)
                    );
                }
                break;
            }
            case OpalSampleMethod::RANDOM_SEQUENCE_UNIFORM: {
                if (Attributes::getReal(itsAttr[SEED])) {
                    sampleMethod_m.reset(
                        new SampleRandomizedSequence<double>(lower, upper, step, seed)
                    );
                } else {
                    sampleMethod_m.reset(
                        new SampleRandomizedSequence<double>(lower, upper, step)
                    );
                }
                break;
            }
            default: {
                throw OpalException("OpalSample::initialize",
                                    "Invalid \"TYPE\" for the \"SAMPLING\" statement");
            }
        }
    }
}


std::string OpalSample::getVariable() const {
    return Attributes::getString(itsAttr[VARIABLE]);
}
