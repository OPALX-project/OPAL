#ifndef classic_src_AbsBeamline_TrackingTest_Tracking_h
#define classic_src_AbsBeamline_TrackingTest_Tracking_h

class Component;

class Tracking
{
public:
    enum var{z, t}; 
    Tracking() {;}
    ~Tracking() {;}

    void integrate(var indie, double target, double* x_in, Component* field, double step_size);

    //Accessors etc
    double getAbsoluteError() const {return _absoluteError;}
    double getRelativeError() const {return _relativeError;}
    int getMaximumStepNumber() const {return _maxNSteps;}
    double getCharge() const {return _q;}
    void setAbsoluteError(double err) {_absoluteError = err;}
    void setRelativeError(double err) {_relativeError = err;}
    void setMaximumStepNumber(int n) {_maxNSteps     = n;}
    void setCharge(double charge) {_q = charge;}

private:
    double _absoluteError;
    double _relativeError;
    int _maxNSteps;
    Component* field_m;
    double _q;
    double _m;

    //single particle equations of motion
    static int t_equations_motion(double t, const double x[8], double dxdt[8], void* params);
    static int z_equations_motion(double z, const double x[8], double dxdt[8], void* params);

    static const double _c_l;
    static Tracking* tracking;
};

#endif
