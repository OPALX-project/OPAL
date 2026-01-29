#ifndef MESHER_H
#define MESHER_H

#include "Algorithms/Vektor.h"
#include "Utilities/MSLang/Triangle.h"

class Mesher {
public:
    Mesher(std::vector<Vector_t> &vertices);

    std::vector<mslang::Triangle> getTriangles() const;

private:
    void orientVerticesCCW();
    bool isConvex(unsigned int i) const;
    double crossProduct(const Vector_t &a,
                        const Vector_t &b) const;
    bool isPointOnLine(unsigned int i,
                       unsigned int j,
                       const Vector_t &pt) const;
    bool isPointRightOfLine(unsigned int i,
                            unsigned int j,
                            const Vector_t &pt) const;
    bool lineSegmentTouchesOrCrossesLine(unsigned int i,
                                         unsigned int j,
                                         unsigned int k,
                                         unsigned int l) const;
    bool isPotentialEdgeIntersected(unsigned int i) const;
    double dotProduct(unsigned int i,
                      unsigned int j,
                      const Vector_t &pt) const;
    bool isPointInsideCone(unsigned int i,
                           unsigned int j,
                           unsigned int jPlusOne,
                           unsigned int jMinusOne) const;
    bool isEar(unsigned int i) const;
    std::vector<unsigned int> findAllEars() const;
    double computeMinimumAngle(unsigned int i) const;
    unsigned int selectBestEar(std::vector<unsigned int> &ears) const;
    void apply();

    std::vector<mslang::Triangle> triangles_m;
    std::vector<Vector_t> vertices_m;
};

#endif