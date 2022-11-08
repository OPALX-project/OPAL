#ifndef CLASSIC_FIELDMAP1DPROFILE2_HH
#define CLASSIC_FIELDMAP1DPROFILE2_HH

#include "Fields/Fieldmap.h"

class _FM1DProfile2: public _Fieldmap {

public:
    virtual ~_FM1DProfile2();

    virtual bool getFieldstrength(const Vector_t &X, Vector_t &strength, Vector_t &info) const;
    virtual bool getFieldDerivative(const Vector_t &X, Vector_t &E, Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);
    virtual void setExitFaceSlope(const double &);
    virtual void setEdgeConstants(const double &bendAngle, const double &entranceAngle, const double &exitAngle);

private:
    _FM1DProfile2(const std::string& filename);

    static FM1DProfile2 create(const std::string& filename);

    virtual void readMap();
    virtual void freeMap();

    double *EngeCoefs_entry_m;
    double *EngeCoefs_exit_m;

    double zbegin_entry_m;
    double zend_entry_m;
    double polynomialOrigin_entry_m;
    int polynomialOrder_entry_m;

    double exit_slope_m;
    double zbegin_exit_m;
    double zend_exit_m;
    double polynomialOrigin_exit_m;
    int polynomialOrder_exit_m;

    bool rectangular_m;

    double length_m;
    double gapHeight_m;

    /// x position in local coordinate system where central trajectory intersects
    /// the exit edge.
    double xExit_m;

    /// z position in local coordinate system where central trajectory intersects
    /// the exit edge.
    double zExit_m;

    /// Cos and sin of the exit edge rotation with respect to the local coordinates.
    double cosExitRotation_m;
    double sinExitRotation_m;

    friend class _Fieldmap;
};

namespace QRDecomposition {
    void solve(double *Matrix, double *Solution, double *rightHandSide, const int &M, const int &N);
}

using FM1DProfile2 = std::shared_ptr<_FM1DProfile2>;

#endif