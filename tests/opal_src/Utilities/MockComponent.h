/*
 *  Copyright (c) 2014, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef UNIT_TESTS_OPAL_SRC_UTILITIES_MOCKCOMPONENT
#define UNIT_TESTS_OPAL_SRC_UTILITIES_MOCKCOMPONENT

#include "BeamlineGeometry/Euclid3DGeometry.h"
#include "AbsBeamline/Component.h"

#include "Algorithms/PartBunchBase.h"

/** Mockups for an Opal Component (e.g. field object). The idea is to test
 *  field lookup routines and placement routines and the like by generating a
 *  "fake" component.
 */
class MockComponent : public Component {
public:
    MockComponent() : Component("MockComponent"), geom_m(NULL) {
        // std::cout << "MOCK CONSTRUCTOR " << this << std::endl;
    }
    MockComponent(const MockComponent& rhs)
        : Component("MockComponent"), geom_m(rhs.geom_m) {
        // std::cout << "MOCK COPY CONSTRUCTOR " << this << std::endl;
    }
    ~MockComponent() { };//std::cout << "MOCK DESTRUCTOR " << this << std::endl;}
    void accept(BeamlineVisitor&) const {}
    ElementBase* clone() const {return new MockComponent(*this);}
    EMField& getField() {EMField* em = NULL; return *em;}
    EMField& getField() const {EMField* em = NULL; return *em;}
    bool apply(const double&, Vector_t&, Vector_t&) {
        return false;
    }
    bool apply(const size_t&, const double&, Vector_t&, Vector_t&) {
        return true;
    }
    bool apply(const Vector_t& r, const Vector_t& /*P*/, const double& /*t*/,
               Vector_t& E, Vector_t& B) {
        if (r(0) < 0. || r(0) > 1. ||
            r(1) < -1. || r(1) > 0. ||
            r(2) < 0. || r(2) > 1.)
            return true; // isOutOfBounds
        B(0) = r(0);
        B(1) = r(1);
        B(2) = r(2);
        E(0) = -r(0);
        E(1) = -r(1);
        E(2) = -r(2);
        return false; // NOT isOutOfBounds
    }
    void initialise(PartBunchBase<double, 3>*, double&, double&) {}
    void finalise() {}
    bool bends() const {return true;}
    void getDimensions(double&, double&) const {}

    Euclid3DGeometry& getGeometry() {return *geom_m;}
    const Euclid3DGeometry& getGeometry() const  {return *geom_m;}

    // caller has responsibility for this memory
    Euclid3DGeometry* geom_m;
private:
};


class MockComponent2 : public Component {
public:
    MockComponent2() : Component("MockComponent"), geom_m(NULL), refB(1,2,3) {
    }
    MockComponent2(const MockComponent2& rhs)
        : Component("MockComponent"), geom_m(rhs.geom_m), refB(rhs.refB) {
    }
    ~MockComponent2() { };
    void accept(BeamlineVisitor&) const {}
    ElementBase* clone() const {return new MockComponent2(*this);}
    EMField& getField() {EMField* em = NULL; return *em;}
    EMField& getField() const {EMField* em = NULL; return *em;}
    bool apply(const double&, Vector_t&, Vector_t&) {
        return false;
    }
    bool apply(const size_t&, const double&, Vector_t&, Vector_t&) {
        return true;
    }
    bool apply(const Vector_t& r, const Vector_t& /*P*/, const double& /*t*/,
               Vector_t& E, Vector_t& B) {
        lastPos = r;
        if (r(0) < -1. || r(0) > 1. ||
            r(1) < -1. || r(1) > 1. ||
            r(2) < 0. || r(2) > 1.)
            return true; // isOutOfBounds
        B = refB;
        E(0) = 0.;
        E(1) = 0.;
        E(2) = 0.;
        return false; // NOT isOutOfBounds
    }
    void initialise(PartBunchBase<double, 3>*, double&, double&) {}
    void finalise() {}
    bool bends() const {return true;}
    void getDimensions(double&, double&) const {}

    Euclid3DGeometry& getGeometry() {return *geom_m;}
    const Euclid3DGeometry& getGeometry() const  {return *geom_m;}

    // caller has responsibility for this memory
    Euclid3DGeometry* geom_m;
    Vector_t refB;
    Vector_t lastPos;
private:
};


#endif