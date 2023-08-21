#ifndef COORDINATESYSTEMTRAFO
#define COORDINATESYSTEMTRAFO

#include "Algorithms/Vektor.h"
#include "Algorithms/Quaternion.h"
#include "Algorithms/Matrix.h"
#include <boost/numeric/ublas/matrix.hpp>

class CoordinateSystemTrafo {
public:
    CoordinateSystemTrafo();

    CoordinateSystemTrafo(const CoordinateSystemTrafo &right);

    CoordinateSystemTrafo(const Vector_t &origin,
                          const Quaternion &orientation);

    Vector_t transformTo(const Vector_t &r) const;
    Vector_t transformFrom(const Vector_t &r) const;

    Vector_t rotateTo(const Vector_t &r) const;
    Vector_t rotateFrom(const Vector_t &r) const;

    void invert();
    CoordinateSystemTrafo inverted() const;

    CoordinateSystemTrafo& operator=(const CoordinateSystemTrafo& right) = default;
    CoordinateSystemTrafo operator*(const CoordinateSystemTrafo &right) const;
    void operator*=(const CoordinateSystemTrafo &right);

    Vector_t getOrigin() const;
    Quaternion getRotation() const;

    void print(std::ostream&) const;
private:
    Vector_t origin_m;
    Quaternion orientation_m;
    matrix_t rotationMatrix_m;
};

inline
std::ostream& operator<<(std::ostream& os, const CoordinateSystemTrafo &trafo) {
    trafo.print(os);
    return os;
}

inline
Inform& operator<<(Inform& os, const CoordinateSystemTrafo &trafo) {
    trafo.print(os.getStream());
    return os;
}

inline
void CoordinateSystemTrafo::print(std::ostream &os) const {
    os << "Origin: " << origin_m << "\n"
       << "z-axis: " << orientation_m.conjugate().rotate(Vector_t(0,0,1)) << "\n"
       << "x-axis: " << orientation_m.conjugate().rotate(Vector_t(1,0,0));
}

inline
Vector_t CoordinateSystemTrafo::getOrigin() const {
    return origin_m;
}

inline
Quaternion CoordinateSystemTrafo::getRotation() const {
    return orientation_m;
}

inline
CoordinateSystemTrafo CoordinateSystemTrafo::inverted() const {
    CoordinateSystemTrafo result(*this);
    result.invert();

    return result;
}

inline Vector_t CoordinateSystemTrafo::transformTo(const Vector_t& r) const {
    boost::numeric::ublas::vector<double> result(3);
    for (size_t i = 0; i < 3; ++i) {
        result(i) = r[i]-origin_m[i];
    }
    result = boost::numeric::ublas::prod(rotationMatrix_m, result);
    Vector_t transformedVector(result(0), result(1), result(2));
    return transformedVector;
}

inline Vector_t CoordinateSystemTrafo::transformFrom(const Vector_t& r) const {
    return rotateFrom(r) + origin_m;
}

inline Vector_t CoordinateSystemTrafo::rotateTo(const Vector_t& r) const {
    boost::numeric::ublas::vector<double> result(3);
    for (size_t i = 0; i < 3; ++i) {
        result(i) = r[i];
    }
    result = boost::numeric::ublas::prod(rotationMatrix_m, result);
    Vector_t rotatedVector(result(0), result(1), result(2));

    return rotatedVector;
}

inline Vector_t CoordinateSystemTrafo::rotateFrom(const Vector_t& r) const {
    boost::numeric::ublas::vector<double> result(3);
    for (size_t i = 0; i < 3; ++i) {
        result(i) = r[i];
    }
    result = boost::numeric::ublas::prod(boost::numeric::ublas::trans(rotationMatrix_m), result);
    Vector_t rotatedVector(result(0), result(1), result(2));
    return rotatedVector;
}

#endif
