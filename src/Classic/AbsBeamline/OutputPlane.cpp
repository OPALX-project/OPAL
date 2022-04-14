#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "AbsBeamline/OutputPlane.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
//#include "Solvers/CollimatorPhysics.hh"



extern Inform *gmsg;
gsl_rng* OutputPlane::rng_m = gsl_rng_alloc(gsl_rng_mt19937);

OutputPlane::OutputPlane() : OutputPlane("")
{}

OutputPlane::OutputPlane(const std::string &name):
    PluginElement(name),
    maxIterations_m(10),
    tolerance_m(1e-9),
    cLight_m(Physics::c*1.0e-9), // metres per ns
    chargeToMassUnits_m(0),
    positionUnits_m(1e3),
    algorithm_m(INTERPOLATION),
    verbose_m(0) {
}

OutputPlane::OutputPlane(const OutputPlane &right):
    PluginElement(right),
    maxIterations_m(10),
    tolerance_m(1e-9),
    cLight_m(Physics::c*1.0e-9), // metres per ns
    chargeToMassUnits_m(0),
    positionUnits_m(1e3) {

    field_m = right.field_m;
    normal_m = right.normal_m;
    centre_m = right.centre_m;
    radialExtent_m = right.radialExtent_m;
    verticalExtent_m = right.verticalExtent_m;
    horizontalExtent_m = right.horizontalExtent_m;
    maxIterations_m = right.maxIterations_m;
    tolerance_m = right.tolerance_m;
    nullfield_m = right.nullfield_m;
    geom_m = right.geom_m;
    recentre_m = right.recentre_m;
    algorithm_m = right.algorithm_m;
    verbose_m = right.verbose_m;
    material_m = right.material_m;
    materialThickness_m = right.materialThickness_m;
}

OutputPlane::~OutputPlane()
{}

ElementBase* OutputPlane::clone() const {
    ElementBase* element = dynamic_cast<ElementBase*>(new OutputPlane(*this));
    return element;
    //return new OutputPlane(*this);
}


void OutputPlane::accept(BeamlineVisitor & visitor) const {
    visitor.visitOutputPlane(*this);
}

void OutputPlane::doInitialise(PartBunchBase<double, 3>* /*bunch*/) {
    //void OutputPlane::doInitialise() {
    *gmsg << "* Initialize OutputPlane at " << centre_m << " with normal " << normal_m;
    if (radialExtent_m) {
        *gmsg << " radial extent " << radialExtent_m;
    }
    if (horizontalExtent_m) {
        *gmsg << " horizontal extent " << horizontalExtent_m;
    }
    if (verticalExtent_m) {
        *gmsg << " vertical extent " << verticalExtent_m;
    }
    *gmsg << endl;
}

void OutputPlane::doGoOffline() {
    *gmsg << "* OutputPlane goes offline " << getName() << endl;
}

void OutputPlane::RK4Step(
                     const double& tStep,
                     const double& chargeToMass,
                     const double& t,
                     Vector_t& R1,
                     Vector_t& P1) const {
    double thalfStep = tStep/2.;
    double tPlusHalf = t+thalfStep;
    double tPlusStep = t+tStep;
    // f = dy/dt
    // f1 = f(x,t)
    double  deriv1[6];

    getDerivatives(R1, P1, t, chargeToMass, deriv1);
    // f2 = f(x+dt*f1/2, t+dt/2 ).
    double  deriv2[6];
    Vector_t R2(R1[0]+thalfStep*deriv1[0], R1[1]+thalfStep*deriv1[1], R1[2]+thalfStep*deriv1[2]);
    Vector_t P2(P1[0]+thalfStep*deriv1[3], P1[1]+thalfStep*deriv1[4], P1[2]+thalfStep*deriv1[5]);
    getDerivatives(R2, P2, tPlusHalf, chargeToMass, deriv2);

    // f3 = f( x+dt*f2/2, t+dt/2 )
    double  deriv3[6];
    Vector_t R3(R1[0]+thalfStep*deriv2[0], R1[1]+thalfStep*deriv2[1], R1[2]+thalfStep*deriv2[2]);
    Vector_t P3(P1[0]+thalfStep*deriv2[3], P1[1]+thalfStep*deriv2[4], P1[2]+thalfStep*deriv2[5]);
    getDerivatives(R3, P3, tPlusHalf, chargeToMass, deriv3);

    // f4 = f(x+dt*f3, t+dt ).
    double  deriv4[6];
    Vector_t R4(R1[0]+tStep*deriv3[0], R1[1]+tStep*deriv3[1], R1[2]+tStep*deriv3[2]);
    Vector_t P4(P1[0]+tStep*deriv3[3], P1[1]+tStep*deriv3[4], P1[2]+tStep*deriv3[5]);
    getDerivatives(R4, P4, tPlusStep, chargeToMass, deriv4);

    // Return x(t+dt) computed from fourth-order R-K.
    R1[0] += (deriv1[0] + deriv4[0] + 2.*(deriv2[0] + deriv3[0]))*tStep/6.;
    R1[1] += (deriv1[1] + deriv4[1] + 2.*(deriv2[1] + deriv3[1]))*tStep/6.;
    R1[2] += (deriv1[2] + deriv4[2] + 2.*(deriv2[2] + deriv3[2]))*tStep/6.;
    P1[0] += (deriv1[3] + deriv4[3] + 2.*(deriv2[3] + deriv3[3]))*tStep/6.;
    P1[1] += (deriv1[4] + deriv4[4] + 2.*(deriv2[4] + deriv3[4]))*tStep/6.;
    P1[2] += (deriv1[5] + deriv4[5] + 2.*(deriv2[5] + deriv3[5]))*tStep/6.;


}

void OutputPlane::getDerivatives(const Vector_t& R,
                            const Vector_t& P,
                            const double& t,
                            const double& chargeToMass,
                            double* yp) const {
    Vector_t externalB, externalE;

    field_m->apply(R, P, t, externalE, externalB);
	//externalB *= 100.;

    double gamma = sqrt(1 + (P[0]*P[0]+P[1]*P[1]+P[2]*P[2]));
    double gammax = (sqrt(P[0]*P[0] +Physics::m_p*Physics::m_p) - Physics::m_p)/(Physics::m_p) +1;
    double gammay = (sqrt(P[1]*P[1] +Physics::m_p*Physics::m_p) - Physics::m_p)/(Physics::m_p) +1;
    double gammaz = (sqrt(P[2]*P[2] +Physics::m_p*Physics::m_p) - Physics::m_p)/(Physics::m_p) +1;
    
    double betax = sqrt(1 - 1/(gammax*gammax));
    double betay = sqrt(1 - 1/(gammay*gammay));
    double betaz = sqrt(1 - 1/(gammaz*gammaz));
    Vector_t beta = P/gamma;

    //std::cout << "beta[0] = "<< beta[0]<<"\n";
    //std::cout << "bet[1] = "<< beta[1]<<"\n";
    //std::cout << "betax = "<< betax<<"\n";
    //std::cout << "betay = "<< betay<<"\n";
    yp[0] = cLight_m*betax; // [m/ns]
    yp[1] = cLight_m*betay; // [m/ns]
    yp[2] = cLight_m*betaz; // [m/ns]
    double kgms2GeVc = (1e9)*Physics::PCHARGE/Physics::c;
    /*
    yp[3] = chargeToMass*(externalE(0) / cLight_m  + 
                         (externalB(2)*beta[1] - externalB(1)*beta[2])); // [1/ns]
    yp[4] = chargeToMass*(externalE(1) / cLight_m  + 
                         (externalB(0)*beta[2] - externalB(2)*beta[0])); // [1/ns]
    yp[5] = chargeToMass*(externalE(2) / cLight_m  + 
                         (externalB(1)*beta[0] - externalB(0)*beta[1])); // [1/ns]
    */

    std::cout << "Ex = " << externalE(0)<<"\n";
    yp[3] = (chargeToMass*(externalE(0)  + (externalB(2)*betay*cLight_m - externalB(1)*betaz*cLight_m)))/kgms2GeVc; // [1/ns]
    yp[4] = (chargeToMass*(externalE(1)  + (externalB(0)*betaz*cLight_m - externalB(2)*betax*cLight_m)))/kgms2GeVc; // [1/ns]
    yp[5] = (chargeToMass*(externalE(2)  + (externalB(1)*betax*cLight_m - externalB(0)*betay*cLight_m)))/kgms2GeVc; // [1/ns]


}

double getHorizontalDistance(Vector_t R, Vector_t centre_m, Vector_t normal_m) {
    double distance = dot(normal_m, R-centre_m);
    Vector_t intersection = R+distance*normal_m; // position of the projection of R along the normal onto the plane
    Vector_t planarDisplacement = intersection-centre_m; // distance between intersection and the centre
    double horizontalDistance = euclidean_norm(planarDisplacement);
    return horizontalDistance;    
}

bool OutputPlane::checkOne(const int index, const double tstep, double chargeToMass, 
                           double& t, Vector_t& R, Vector_t& P) {

    // distance from particle to the output plane
    Vector_t delta = R-centre_m;
    double distance = dot(normal_m, delta);
    double sStep = tstep*euclidean_norm(P)*cLight_m;
    output_m.str(""); 
    output_m << "* Check crossing of plane " << getName() << " at " << centre_m << " with normal " << normal_m << std::endl;
    output_m << "    Particle " << index << " with R " << R << " P " << P << " tstep " << tstep << std::endl;
    output_m << "    Distance prestep " << distance << " compared to s step length at c " << sStep << std::endl;
    if (fabs(distance) > sStep) {
        // we can't cross the plane
        return false;
    }

    /*
    double horizontalDistance = getHorizontalDistance(R, centre_m, normal_m);
    output_m << "    Horizontal extent " << horizontalExtent_m << " distance " << horizontalDistance
             << "  compared to max horizontal step " << (horizontalExtent_m+sStep)*(horizontalExtent_m+sStep) << std::endl;
    if (horizontalDistance >
        (horizontalExtent_m+sStep)*(horizontalExtent_m+sStep)) {
        // we can't make the aperture even travelling along the plane at c_light
        return false;
    }*/
    Vector_t rTest(R);
    Vector_t pTest(P);
    RK4Step(tstep, chargeToMass, t, rTest, pTest);
    double distanceTest = dot(normal_m, (rTest-centre_m));
    output_m << "    Test step pos: " << rTest << " P " << pTest
             << " distance prestep " << distance << " distance postStep " << distanceTest << std::endl;

    if (distanceTest/distance > 0) {
        // step does not cross the plane - give up
        // note that particle could cross and cross back the plane in a single
        // time-step; in this case the particle is not registered
        return false;
    }

    if (algorithm_m == RK4) {

        rk4Test(tstep, chargeToMass, t, R, P);
    } else if (algorithm_m == INTERPOLATION) {
        interpolation(t, R, P);
    }
    delta = R-centre_m;
    //horizontalDistance = getHorizontalDistance(R, centre_m, normal_m);
    //if (horizontalExtent_m > 0 && horizontalDistance >
    //    (horizontalExtent_m+sStep)*(horizontalExtent_m+sStep)) {
    //    // out of horizontal extent
    //    return false;
    // }
    if (horizontalExtent_m > 0 && delta[0]*delta[0]+delta[1]*delta[1] > horizontalExtent_m*horizontalExtent_m) {
        // out of horizontal extent. Again, defined in global coordinates
        return false;
    }
    if (verticalExtent_m > 0 && abs(delta[2]) > verticalExtent_m) {
        // out of vertical extent; note that this is defined in global coordinates
        // not in the coordinates of the plane
        return false;
    }
    if (recentre_m == index) {
        recentre(R, P);
        *gmsg << "* Recentred output plane to " << centre_m
            << " with normal " << normal_m << " by event " << index << endl;
    }
    return true;
}

void OutputPlane::rk4Test(double tstep, double chargeToMass,
                          double& t, Vector_t& R, Vector_t& P) {
    double preStepDistance = 0.0;
    double postStepDistance = dot(normal_m, (R-centre_m));
    size_t iteration = 0;
    while (fabs(postStepDistance) > tolerance_m && iteration < maxIterations_m) {
        preStepDistance = postStepDistance;
        RK4Step(tstep, chargeToMass, t, R, P); // this updates R and P
        t += tstep;
        postStepDistance = dot(normal_m, (R-centre_m));
        if (postStepDistance/preStepDistance < 0) { // straddling the plane
            // we stepped too far; step in opposite direction
            // step length in ratio of distance to plane
            tstep *= -abs(postStepDistance)/abs(postStepDistance - preStepDistance);
        } else {
            // we didn't step far enough; step in same direction
            tstep *= abs(postStepDistance)/abs(postStepDistance - preStepDistance);
        }
        output_m << "    Output " << postStepDistance << " R " << R << " P " << P 
                 << " centre " << centre_m << " d: " << R-centre_m << std::endl;
        iteration++;
    }
}

void OutputPlane::interpolation(double& t, Vector_t& R, Vector_t& P) {
    // trajectory R = R0 + V dt 
    // plane = (X-X0).n
    // intersection time t0 = (X0-R0).N / (V.N)
    // intersection position = R0 + V dt
    // relativistic gamma:
    output_m << "    Interpolated from R " << R;
    double gamma = ::sqrt(1+P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
    Vector_t velocity = P/gamma*cLight_m; // m/ns
    double dt = dot((centre_m-R), normal_m) / dot(velocity, normal_m);
    R += velocity*dt;
    t += dt;
    output_m << " to R " << R << " with P " << P << std::endl;
}


bool OutputPlane::doCheck(PartBunchBase<double, 3> *bunch, const int /*turnnumber*/, 
                          const double t, const double tstep) {
    size_t tempnum = bunch->getLocalNum();
    for(unsigned int i = 0; i < tempnum; ++i) {
        Vector_t R(bunch->R[i]/positionUnits_m);
        Vector_t P(bunch->P[i]);
        double t0(t);
        double chargeToMass = bunch->Q[i]/bunch->M[i]*cLight_m*cLight_m/Physics::q_e/1e3;
        bool crossing = checkOne(i, tstep, chargeToMass, t0, R, P);
        if (verbose_m > 9 && (crossing || verbose_m > 19)) {
            *gmsg << output_m.str();
        }
        if (crossing && lossDs_m) {
            nHits_m += 1;
            //lossDs_m->addParticle(R, P, bunch->ID[i], t0, turnnumber);
        }
        if (crossing && materialThickness_m > 0.0) {
            hardCodeMaterialPhysicsGraphite(bunch->P[i]);
        }
    }
    return false;
}
/*
// never used; ccollimator was giving spurious answers
void OutputPlane::applyMaterialPhysics(Vector_t& R, Vector_t& P) const {
    // we set the "effective time step" as thickness/beta/c where beta is 
    // calculated in the direction perpendicular to the foil
    std::cerr << "APPLYMATEIRAL PHYSICS" << std::endl;
    double ns2s = 1e-9;
    double gamma = ::sqrt(1+P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
    double beta = abs(P[0]*normal_m[0]+P[1]*normal_m[1]+P[2]*normal_m[2])/gamma;
    double mass = Physics::m_p;
    double dt = abs(materialThickness_m/beta/cLight_m*ns2s);
    material_m->computeCoulombScattering(R, P, dt);
    double newKineticEnergy = (gamma-1)*mass;
    std::cerr << "APPLY MATERIAL PHYSICS " << (gamma-1)*mass << std::endl;
    // BUG: fluctuations are set to false because it is returning NaN.
    material_m->computeEnergyLoss(newKineticEnergy, dt, false);
    if (newKineticEnergy <= 0) {
        P = Vector_t(0, 0, 0);
    } else {
        double new_gamma = newKineticEnergy/mass+1;
        double new_bg = ::sqrt(new_gamma*new_gamma-1);
        P *= new_bg/beta/gamma;
    }
    std::cerr << "    dt " << dt << " beta " << beta << " gamma " << gamma << " mass " << mass << " new KE " << newKineticEnergy << std::endl;
}
*/

void OutputPlane::recentre(Vector_t R, Vector_t P) {
    setCentre(R);
    setNormal(P);
    recentre_m = -1; // don't recentre again
}
/*
ElementBase::ElementType OutputPlane::getType() const {
    return OUTPUTPLANE;
}
*/
ElementType OutputPlane::getType() const {
    return ElementType::OUTPUTPLANE;
}


/*
  "carbon":{
    "radiation_length": 42.70 ,
    "density": 2.000 ,
    "ion_energy": 78.0,
    "molecular_mass":12.0107,
    "delta":0.10,
    "sigma_0,1":{"E_t":1.36e-5, "a_1":6.54e1, "a_2":1.53, "a_3":1.98,
                                "a_4":-1.8e-2, "a_5":1.15e1, "a_6":1.08},
    "sigma_-1,0":{"E_t":6.4e-7, "a_1":8.58, "a_2":1.32e-1, "a_3":2.61e1,
                                "a_4":8.4e-1, "a_5":1.5e3, "a_6":2.1},
    "elements_list":[
        {"name":"carbon", "z":6, "a":12.0107, "fraction":1.0}
    ],
    "specific_heat":0.71
  },
*/

Vector_t perpendicular(Vector_t a, Vector_t b) {
    Vector_t perp(
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    );
    perp = perp / ::sqrt(perp[0]*perp[0]+perp[1]*perp[1]+perp[2]*perp[2]);
    return perp;
}

void OutputPlane::hardCodeMaterialPhysicsGraphite(Vector_t& P) const {
    double density = 2.0; // g/cm^3
    double ion_energy = 78.0*1e-6; // MeV
    //double delta = 0.10;
    //double molecular_mass = 12.0107; //
    double X0 = 42.70/density; // cm
    double z = 6.0;
    double a = 12.0107;

    double charge = 1.0;
    double gamma = ::sqrt(1+P[0]*P[0]+P[1]*P[1]+P[2]*P[2]);
    double beta = ::sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2])/gamma;
    double mass = Physics::m_p*1e3;
    double m_e = Physics::m_e*1e3;
    double ds = materialThickness_m*100.0; // [cm] assume perpendicular

    // bethe model for beta gamma > 0.05 (proton energy > 1 MeV)
    double bg2 = beta*beta*gamma*gamma;
    double m_ratio = m_e/mass;
    double t_max = 2.0*m_e*bg2/(1.0 + 2.0*gamma*m_ratio+m_ratio*m_ratio);
    double coefficient = log(2*m_e*bg2*t_max/(ion_energy*ion_energy))/2.;
    coefficient -= beta*beta;
    coefficient *= 0.307075*charge*charge/beta/beta*density;

    double dE = coefficient * z/a * ds;
    double gamma_new = gamma - dE/mass;

    double t0 = 13.6/beta/mass/beta/gamma;
    t0 *= abs(charge)*::sqrt(ds/X0);

    // incident particle has some arbitrary direction
    // outgoing particle gets two scatters along vectors perpendicular to the
    // direction of travel and perpendicular to each other. First vector v1 is
    // defined as in perpendicular to z axis zhat and direction of travel P. Second
    // vector v2 is perpendicular to v1 and p
    double xp = gsl_ran_gaussian(rng_m, t0);
    double yp = gsl_ran_gaussian(rng_m, t0);
    Vector_t zhat(0, 0, 1);
    Vector_t v1 = perpendicular(P, zhat)*xp*beta*gamma;
    Vector_t v2 = perpendicular(P, v1)*yp*beta*gamma;
    P += v1;
    P += v2;
    //std::cerr << "hardCodeMaterialPhysicsGraphite " << xp << " " << yp << " " << P << " " << v1 << " " << v2 << std::endl;
    P *= ::sqrt(gamma_new*gamma_new-1)/::sqrt(P[0]*P[0]+P[1]*P[1]+P[2]*P[2]); // beta gamma = (gamma^2-1)^0.5
}

