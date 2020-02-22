#ifndef CLASSIC_GENERALFIELDMAP_H
#define CLASSIC_GENERALFIELDMAP_H

#include "Fields/Fieldmap.h"

#include <gsl/gsl_interp.h>

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#include "tests/classic_src/Fields/GeneralFieldMapTest.cpp"
#endif

class GeneralFieldMap: public Fieldmap {

public:
    virtual bool getFieldstrength(const Vector_t &R, ComplexVector_t &E, ComplexVector_t &B) const;
    virtual bool getFieldDerivative(const Vector_t &R, Vector_t &E,
                                    Vector_t &B, const DiffDirection &dir) const;
    virtual void getFieldDimensions(double &zBegin, double &zEnd) const;
    virtual void getFieldDimensions(double &xIni, double &xFinal,
                                    double &yIni, double &yFinal,
                                    double &zIni, double &zFinal) const;
    virtual void swap();
    virtual void getInfo(Inform *);
    virtual double getFrequency() const;
    virtual void setFrequency(double freq);

    virtual bool isInside(const Vector_t &r) const;

private:
#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GeneralFieldMapTest, ReadInFromElement);
#endif

    GeneralFieldMap(std::string const& aFilename);
    ~GeneralFieldMap();

    virtual void readMap();
    struct FieldData
    {
        double accuracy;
        double frequency;
        std::vector<double> zValues;
        std::vector<double> EzRealValues;
        std::vector<double> EzImagValues;
    };
    FieldData readFieldData() const;
    virtual void freeMap();
    void allocateMemory(unsigned int numGridPoints);
    void evenlyDistributeData(std::vector<double>& zValues,
                              std::vector<double>& EzRealValues,
                              std::vector<double>& EzImagValues) const;
    void initInterpolation(std::vector<double> const& zValues,
                           std::pair<std::vector<double>, std::vector<double>> const& fourCoefs);

    void computeFieldOnAxis(double z,
                            std::vector<std::pair<double, double>> &fieldComponents) const;

    void computeFieldOffAxis(const Vector_t &R,
                             ComplexVector_t &E,
                             ComplexVector_t &B,
                             std::vector<std::pair<double, double>> const& fieldComponents) const;


    double frequency_m;                     /// Field angular frequency (Hz).
    double twoPiOverLambdaSq_m;             /// 2 Pi divided by the field RF wavelength squared.

    double zBegin_m;                        /// Longitudinal start of field.
    double length_m;                        /// Field length.
    unsigned int numberOfGridPoints_m;      /// Number of grid points in field input file.
    double deltaZ_m;                        /// Field grid point spacing.
    unsigned int accuracy_m;

    struct SplineFree {
        void operator() (gsl_spline *ptr) {
            gsl_spline_free(ptr);
        }
    };

    typedef std::unique_ptr<gsl_spline, SplineFree> splinePtr_t;
    class complexSpline: public std::pair<splinePtr_t, splinePtr_t> {
    public:
        splinePtr_t& real;
        splinePtr_t& imag;
        complexSpline():
            std::pair<splinePtr_t, splinePtr_t>(),
            real(first),
            imag(second)
        { }
    };

    std::vector<std::complex<double>> onAxisField_m; /// On axis field data.
    complexSpline onAxisFieldInterpolants_m;       /// On axis field interpolation structure.
    complexSpline onAxisFieldPInterpolants_m;      /// On axis field first derivative interpolation structure.
    complexSpline onAxisFieldPPInterpolants_m;     /// On axis field second derivative interpolation structure.
    complexSpline onAxisFieldPPPInterpolants_m;    /// On axis field third derivative interpolation structure.

    struct InterpAccelFree {
        void operator() (gsl_interp_accel *ptr) {
            gsl_interp_accel_free(ptr);
        }
    };
    /// Corresponding interpolation evaluation accelerators.
    std::unique_ptr<gsl_interp_accel, InterpAccelFree> onAxisFieldAccel_m;
    // std::unique_ptr<gsl_interp_accel, gsl_interp_accel_free> onAxisFieldPAccel_m;
    // std::unique_ptr<gsl_interp_accel, gsl_interp_accel_free> onAxisFieldPPAccel_m;
    // std::unique_ptr<gsl_interp_accel, gsl_interp_accel_free> onAxisFieldPPPAccel_m;

    friend class Fieldmap;
};

inline bool GeneralFieldMap::isInside(const Vector_t &r) const
{
    return r(2) >= zBegin_m && r(2) < zBegin_m + length_m;
}

#endif