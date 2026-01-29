#include "Utilities/Mesher.h"

Mesher::Mesher(std::vector<Vector_t> &vertices):
    vertices_m(vertices)
{
    apply();
}

std::vector<mslang::Triangle> Mesher::getTriangles() const {
    return triangles_m;
}

void Mesher::orientVerticesCCW() {
    const unsigned int size = vertices_m.size();
    double sum = 0.0;
    for (unsigned int i = 0; i < size; ++ i) {
        unsigned int iPlusOne = (i + 1) % size;
        Vector_t edge(vertices_m[iPlusOne][0] - vertices_m[i][0],
                      vertices_m[iPlusOne][1] + vertices_m[i][1],
                      0.0);
        sum += edge[0] * edge[1];
    }
    if (sum <= 0.0) return;

    std::vector<Vector_t> reoriented(vertices_m.rbegin(), vertices_m.rend());

    vertices_m.swap(reoriented);
}

bool Mesher::isConvex(unsigned int i) const {
    unsigned int numVertices = vertices_m.size();
    i = i % numVertices;
    unsigned int iPlusOne = (i + 1) % numVertices;
    unsigned int iMinusOne = (i + numVertices - 1) % numVertices;

    Vector_t edge0 = vertices_m[iMinusOne] - vertices_m[i];
    Vector_t edge1 = vertices_m[iPlusOne] - vertices_m[i];

    double vectorProduct = edge0[0] * edge1[1] - edge0[1] * edge1[0];

    return (vectorProduct < 0.0);
}

double Mesher::crossProduct(const Vector_t &a,
                                 const Vector_t &b) const {
    return a[0] * b[1] - a[1] * b[0];
}

bool Mesher::isPointOnLine(unsigned int i,
                                unsigned int j,
                                const Vector_t &pt) const {
    Vector_t aTmp = vertices_m[j] - vertices_m[i];
    Vector_t bTmp = pt - vertices_m[i];

    double r = crossProduct(aTmp, bTmp);

    return std::abs(r) < 1e-10;
}

bool Mesher::isPointRightOfLine(unsigned int i,
                                unsigned int j,
                                const Vector_t &pt) const {
    Vector_t aTmp = vertices_m[j] - vertices_m[i];
    Vector_t bTmp = pt - vertices_m[i];

    return crossProduct(aTmp, bTmp) < 0.0;
}

bool Mesher::lineSegmentTouchesOrCrossesLine(unsigned int i,
                                                  unsigned int j,
                                                  unsigned int k,
                                                  unsigned int l) const {
    return (isPointOnLine(i, j, vertices_m[k]) ||
            isPointOnLine(i, j, vertices_m[l]) ||
            (isPointRightOfLine(i, j, vertices_m[k]) ^
             isPointRightOfLine(i, j, vertices_m[l])));
}

bool Mesher::isPotentialEdgeIntersected(unsigned int i) const {
    unsigned int numVertices = vertices_m.size();
    if (numVertices < 5) return false;

    i = i % numVertices;
    unsigned int iPlusOne = (i + 1) % numVertices;
    unsigned int iMinusOne = (i + numVertices - 1) % numVertices;

    mslang::BoundingBox bbPotentialNewEdge;
    bbPotentialNewEdge.center_m = 0.5 * (vertices_m[iMinusOne] + vertices_m[iPlusOne]);
    bbPotentialNewEdge.width_m = std::abs(vertices_m[iMinusOne][0] - vertices_m[iPlusOne][0]);
    bbPotentialNewEdge.height_m = std::abs(vertices_m[iMinusOne][1] - vertices_m[iPlusOne][1]);

    for (unsigned int j = iPlusOne + 1; j < iPlusOne + numVertices - 3; ++ j) {
        unsigned int k = (j % numVertices);
        unsigned int kPlusOne = ((k + 1) % numVertices);

        mslang::BoundingBox bbThisEdge;
        bbThisEdge.center_m = 0.5 * (vertices_m[k] + vertices_m[kPlusOne]);
        bbThisEdge.width_m = std::abs(vertices_m[k][0] - vertices_m[kPlusOne][0]);
        bbThisEdge.height_m = std::abs(vertices_m[k][1] - vertices_m[kPlusOne][1]);

        if (bbPotentialNewEdge.doesIntersect(bbThisEdge) &&
            lineSegmentTouchesOrCrossesLine(iPlusOne, iMinusOne,
                                            k, kPlusOne) &&
            lineSegmentTouchesOrCrossesLine(k, kPlusOne,
                                            iPlusOne, iMinusOne))
            return true;
    }

    return false;
}

double getAngleBetweenEdges(const Vector_t &a,
                            const Vector_t &b) {
    double lengthA = mslang::euclidean_norm2D(a);
    double lengthB = mslang::euclidean_norm2D(b);

    double angle = std::acos((a[0] * b[0] + a[1] * b[1]) / lengthA / lengthB);
    if (a[0] * b[1] - a[1] * b[0] < 0.0)
        angle += M_PI;

    return angle;
}

double Mesher::dotProduct(unsigned int i,
                          unsigned int j,
                          const Vector_t &pt) const {
    Vector_t edge0 = vertices_m[j] - vertices_m[i];
    Vector_t edge1 = vertices_m[j] - pt;

    return edge0[0] * edge1[0] + edge0[1] * edge1[1];
}

double dotProduct(const Vector_t &a,
                  const Vector_t &b) {
    return a[0] * b[0] + a[1] * b[1];
}

bool Mesher::isPointInsideCone(unsigned int i,
                               unsigned int j,
                               unsigned int jPlusOne,
                               unsigned int jMinusOne) const {
    const Vector_t &pt = vertices_m[i];

    return !((isPointRightOfLine(jMinusOne, j, pt) &&
              dotProduct(jMinusOne, j, pt) > 0.0) ||
             (isPointRightOfLine(j, jPlusOne, pt) &&
              dotProduct(jPlusOne, j, pt) < 0.0));
}

bool Mesher::isEar(unsigned int i) const {
    unsigned int size = vertices_m.size();

    unsigned int iMinusTwo = (i + size - 2) % size;
    unsigned int iMinusOne = (i + size - 1) % size;
    unsigned int iPlusOne = (i + 1) % size;
    unsigned int iPlusTwo = (i + 2) % size;

    bool convex = isConvex(i);
    bool isInsideCone1 = isPointInsideCone(iMinusOne, iPlusOne, iPlusTwo, i);
    bool isInsideCone2 = isPointInsideCone(iPlusOne, iMinusOne, i, iMinusTwo);
    bool notCrossed = !isPotentialEdgeIntersected(i);

    return (convex &&
            isInsideCone1 &&
            isInsideCone2 &&
            notCrossed);
}

std::vector<unsigned int> Mesher::findAllEars() const {
    unsigned int size = vertices_m.size();
    std::vector<unsigned int> ears;

    for (unsigned int i = 0; i < size; ++ i) {
        if (isEar(i)) {
            ears.push_back(i);
        }
    }

    return ears;
}

double Mesher::computeMinimumAngle(unsigned int i) const {
    unsigned int numVertices = vertices_m.size();
    unsigned int previous = (i + numVertices - 1) % numVertices;
    unsigned int next = (i + 1) % numVertices;

    Vector_t edge0 = vertices_m[i] - vertices_m[previous];
    Vector_t edge1 = vertices_m[next] - vertices_m[i];
    Vector_t edge2 = vertices_m[previous] - vertices_m[next];
    double length0 = mslang::euclidean_norm2D(edge0);
    double length1 = mslang::euclidean_norm2D(edge1);
    double length2 = mslang::euclidean_norm2D(edge2);

    double angle01 = std::acos(-(edge0[0] * edge1[0] + edge0[1] * edge1[1]) / length0 / length1);
    double angle12 = std::acos(-(edge1[0] * edge2[0] + edge1[1] * edge2[1]) / length1 / length2);
    double angle20 = M_PI - angle01 - angle12;

    return std::min(std::min(angle01, angle12), angle20);
}

unsigned int Mesher::selectBestEar(std::vector<unsigned int> &ears) const {
    unsigned int numEars = ears.size();

    double maxMinAngle = computeMinimumAngle(ears[0]);
    unsigned int earWithMaxMinAngle = 0;

    for (unsigned int i = 1; i < numEars; ++ i) {
        double angle = computeMinimumAngle(ears[i]);

        if (angle > maxMinAngle) {
            maxMinAngle = angle;
            earWithMaxMinAngle = i;
        }
    }

    return ears[earWithMaxMinAngle];
}

void Mesher::apply() {
    orientVerticesCCW();

    unsigned int numVertices = vertices_m.size();
    while (numVertices > 3) {
        std::vector<unsigned int> allEars = findAllEars();
        unsigned int bestEar = selectBestEar(allEars);
        unsigned int next = (bestEar + 1) % numVertices;
        unsigned int previous = (bestEar + numVertices - 1) % numVertices;

        mslang::Triangle tri;
        tri.nodes_m[0] = vertices_m[previous];
        tri.nodes_m[1] = vertices_m[bestEar];
        tri.nodes_m[2] = vertices_m[next];
        // tri.print(4);
        triangles_m.push_back(tri);

        vertices_m.erase(vertices_m.begin() + bestEar);

        -- numVertices;
    }

    mslang::Triangle tri;
    tri.nodes_m = vertices_m;
    triangles_m.push_back(tri);
}