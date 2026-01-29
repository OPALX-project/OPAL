//
// Mesh Generator
//
// Copyright (c) 2008-2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include "Algorithms/Vektor.h"

#include <string>
#include <utility>
#include <vector>

class ElementBase;

class MeshData {
public:
    std::vector<Vector_t> vertices_m;
    std::vector<Vektor<unsigned int, 3> > triangles_m;
    std::vector<std::pair<Vector_t, Vector_t> > decorations_m;
    int type_m;
};


class MeshGenerator {
public:
    MeshGenerator();

    void add(const ElementBase &element);

    void write(const std::string &fname);
private:
    enum MeshType {OTHER = 0
                 , DIPOLE
                 , QUADRUPOLE
                 , SEXTUPOLE
                 , OCTUPOLE
                 , SOLENOID
                 , RFCAVITY
    };

    static MeshData getCylinder(double length,
                                double minor,
                                double major,
                                double formFactor,
                                const unsigned int numSegments = 36);

    static MeshData getBox(double length,
                           double width,
                           double height,
                           double formFactor);


    std::vector<MeshData> elements_m;
};

#endif
