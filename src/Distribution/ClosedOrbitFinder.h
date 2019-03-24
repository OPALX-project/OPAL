/**
 * @file ClosedOrbitFinder.h
 * The algorithm is based on the paper of M. M. Gordon: "Computation of closed orbits and basic focusing properties for
 * sector-focused cyclotrons and the design of 'cyclops'" (1983)
 * As template arguments one chooses the type of the variables and the integrator for the ODEs. The supported steppers can
 * be found on
 * http://www.boost.org/ where it is part of the library Odeint.
 *
 * @author Matthias Frey
 * @version 1.0
 */

#ifndef CLOSEDORBITFINDER_H
#define CLOSEDORBITFINDER_H

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"

#include "AbstractObjects/OpalData.h"

#include "AbsBeamline/Cyclotron.h"

// include headers for integration
#include <boost/numeric/odeint/integrate/integrate_n_steps.hpp>
#include <boost/filesystem.hpp>

extern Inform *gmsg;

/// Finds a closed orbit of a cyclotron for a given energy
template<typename Value_type, typename Size_type, class Stepper>
class ClosedOrbitFinder
{
    public:
        /// Type of variables
        typedef Value_type value_type;
        /// Type for specifying sizes
        typedef Size_type size_type;
        /// Type of container for storing quantities (path length, orbit, etc.)
        typedef std::vector<value_type> container_type;
        /// Type for holding state of ODE values
        typedef std::vector<value_type> state_type;

        typedef std::function<void(const state_type&, state_type&, const double)> function_t;

        /// Sets the initial values for the integration and calls findOrbit().
        /*!
         * @param E0 is the potential energy (particle energy at rest) [MeV].
         * @param N specifies the number of splits (2pi/N), i.e number of integration steps
         * @param cycl is the cyclotron element
         * @param domain is a boolean (default: true). If "true" the closed orbit is computed over a single sector,
         * otherwise over 2*pi.
         * @param Nsectors is an int (default: 1). Number of sectors that the field map is averaged over
         * in order to avoid first harmonic. Only valid if domain is false
         */
        ClosedOrbitFinder(value_type E0, size_type N,
                          Cyclotron* cycl, bool domain = true, int Nsectors = 1);

        /// Returns the inverse bending radius (size of container N+1)
        container_type getInverseBendingRadius(const value_type& angle = 0);

        /// Returns the step lengths of the path (size of container N+1)
        container_type getPathLength(const value_type& angle = 0);

        /// Returns the field index (size of container N+1)
        container_type getFieldIndex(const value_type& angle = 0);

        /// Returns the radial and vertical tunes (in that order)
        std::pair<value_type,value_type> getTunes();

        /// Returns the closed orbit (size of container N+1) starting at specific angle (only makes sense when computing
        /// the closed orbit for a whole turn) (default value: 0°).
        /// Attention: It computes the starting index of the array. If it's not an integer it just cuts the floating point
        /// part, i.e. it takes the next starting index below. There's no interpolation of the radius.
        /*!
         * @param angle is the start angle for the output. Has to be within [0°,360°[ (default: 0°).
         */
        container_type getOrbit(value_type angle = 0);

        /// Returns the momentum of the orbit (size of container N+1)starting at specific angle (only makes sense when
        /// computing the closed orbit for a whole turn) (default value: 0°), \f$ \left[ p_{r} \right] = \si{m}\f$.
        /// Attention: It computes the starting index of the array. If it's not an integer it just cuts the floating point
        /// part, i.e. it takes the next starting index below. There's no interpolation of the momentum.
        /*!
         * @param angle is the start angle for the output. Has to be within [0°,360°[ (default: 0°).
         * @returns the momentum in \f$ \beta * \gamma \f$ units
         */
        container_type getMomentum(value_type angle = 0);

        /// Returns the average orbit radius
        value_type getAverageRadius();

        /// Returns the frequency error
        value_type getFrequencyError();

        /// Returns true if a closed orbit could be found
        bool isConverged();

        /// Computes the closed orbit
        /*!
         * @param accuracy specifies the accuracy of the closed orbit
         * @param maxit is the maximal number of iterations done for finding the closed orbit
         * @param ekin energy for which to find closed orbit (in tune mode: upper limit of range)
         * @param dE step increase [MeV]
         * @param rguess initial radius guess in [mm]
         * @param isTuneMode comptute tunes of all energies in one sweep
         */
        bool findOrbit(value_type accuracy, size_type maxit,
                       value_type ekin,
                       value_type dE = 1.0, value_type rguess = -1.0,
                       bool isTuneMode = false);

        /// Fills in the values of h_m, ds_m, fidx_m.
        void computeOrbitProperties(const value_type& E);

    private:
        /// This function is called by the function getTunes().
        /*! Transfer matrix Y = [y11, y12; y21, y22] (see Gordon paper for more details).
         * @param y are the positions (elements y11 and y12 of Y)
         * @param py2 is the momentum of the second solution (element y22 of Y)
         * @param ncross is the number of sign changes (\#crossings of zero-line)
         */
        value_type computeTune(const std::array<value_type,2>&, value_type, size_type);

        // Compute closed orbit for given energy
        bool findOrbitOfEnergy_m(const value_type&, container_type&, value_type&,
                                 const value_type&, size_type);

        /// This function computes nzcross_ which is used to compute the tune in z-direction and the frequency error
//         void computeVerticalOscillations();

        /// This function rotates the calculated closed orbit finder properties to the initial angle
        container_type rotate(value_type angle, container_type& orbitProperty);

        /// Stores current position in horizontal direction for the solutions of the ODE with different initial values
        std::array<value_type,2> x_m; // x_m = [x1, x2]
        /// Stores current momenta in horizontal direction for the solutions of the ODE with different initial values
        std::array<value_type,2> px_m; // px_m = [px1, px2]
        /// Stores current position in vertical direction for the solutions of the ODE with different initial values
        std::array<value_type,2> z_m; // z_m = [z1, z2]
        /// Stores current momenta in vertical direction for the solutions of the ODE with different initial values
        std::array<value_type,2> pz_m; // pz_m = [pz1, pz2]

        /// Stores the inverse bending radius
        container_type h_m;
        /// Stores the step length
        container_type ds_m;
        /// Stores the radial orbit (size: N_m+1)
        container_type r_m;
        /// Stores the vertical oribt (size: N_m+1)
        container_type vz_m;
        /// Stores the radial momentum
        container_type pr_m;
        /// Stores the vertical momentum
        container_type vpz_m;
        /// Stores the field index
        container_type fidx_m;

        /// Counts the number of zero-line crossings in horizontal direction (used for computing horizontal tune)
        size_type nxcross_m;
        /// Counts the number of zero-line crossings in vertical direction (used for computing vertical tune)
        size_type nzcross_m; //#crossings of zero-line in x- and z-direction

        /// Is the rest mass [MeV / c**2]
        value_type E0_m;

        /// Is the nominal orbital frequency
        /* (see paper of Dr. C. Baumgarten: "Transverse-Longitudinal
         * Coupling by Space Charge in Cyclotrons" (2012), formula (1))
         */
        value_type wo_m;
        /// Number of integration steps
        size_type N_m;
        /// Is the angle step size
        value_type dtheta_m;

        /// Is the average radius
        value_type ravg_m;

        /// Is the phase
        value_type phase_m;

        /**
         * Stores the last orbit value (since we have to return to the beginning to check the convergence in the
         * findOrbit() function. This last value is then deleted from the array but is stored in lastOrbitVal_m to
         * compute the vertical oscillations)
         */
        /* value_type lastOrbitVal_m; */

        /**
         * Stores the last momentum value (since we have to return to the beginning to check the convergence in the
         * findOrbit() function. This last value is then deleted from the array but is stored in lastMomentumVal_m to
         * compute the vertical oscillations)
         */
        /* value_type lastMomentumVal_m; */

        /**
         * Boolean which is true by default. "true": orbit integration over one sector only, "false": integration
         * over 2*pi
         */
        bool domain_m;
        /**
         * Number of sectors to average the field map over
         * in order to avoid first harmonic. Only valid if domain is false
         */
        int  nSectors_m;

        /// Defines the stepper for integration of the ODE's
        Stepper stepper_m;

        /*!
         * This quantity is defined in the paper "Transverse-Longitudinal Coupling by Space Charge in Cyclotrons"
         * of Dr. Christian Baumgarten (2012)
         * The lambda function takes the orbital frequency \f$ \omega_{o} \f$ (also defined in paper) as input argument.
         */
        std::function<double(double)> acon_m = [](double wo) { return Physics::c / wo; };

        /// Cyclotron unit \f$ \left[T\right] \f$ (Tesla)
        /*!
         * The lambda function takes the orbital frequency \f$ \omega_{o} \f$ as input argument.
         */
        std::function<double(double, double)> bcon_m = [](double e0, double wo) {
            return e0 * 1.0e7 / (/* physics::q0 */ 1.0 * Physics::c * Physics::c / wo);
        };

        Cyclotron* cycl_m;
};

// -----------------------------------------------------------------------------------------------------------------------
// PUBLIC MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------

template<typename Value_type, typename Size_type, class Stepper>
ClosedOrbitFinder<Value_type,
                  Size_type,
                  Stepper>::ClosedOrbitFinder(value_type E0,
                                              size_type N, Cyclotron* cycl,
                                              bool domain, int nSectors)
    : nxcross_m(0)
    , nzcross_m(0)
    , E0_m(E0)
    , wo_m(cycl->getRfFrequ()*1E6/cycl->getCyclHarm()*2.0*Physics::pi)
    , N_m(N)
    , dtheta_m(Physics::two_pi/value_type(N))
    , ravg_m(0)
    , phase_m(0)
    /* , lastOrbitVal_m(0.0) */
    /* , lastMomentumVal_m(0.0) */
    , domain_m(domain)
    , nSectors_m(nSectors)
    , stepper_m()
    , cycl_m(cycl)
{

    if ( cycl_m->getFMLowE() > cycl_m->getFMHighE() )
        throw OpalException("ClosedOrbitFinder::ClosedOrbitFinder()",
                            "Incorrect cyclotron energy (MeV) bounds: Maximum cyclotron energy smaller than minimum cyclotron energy.");

    // if domain_m = true --> integrate over a single sector
    if (domain_m) {
        N_m /=  cycl_m->getSymmetry();
    }

    cycl_m->read(cycl_m->getFieldFlag(cycl_m->getCyclotronType()),
                 cycl_m->getBScale());

    // reserve storage for the orbit and momentum (--> size = 0, capacity = N_m+1)
    /*
     * we need N+1 storage, since dtheta = 2pi/N (and not 2pi/(N-1)) that's why we need N+1 integration steps
     * to return to the origin (but the return size is N_m)
     */
    r_m.reserve(N_m + 1);
    pr_m.reserve(N_m + 1);
    vz_m.reserve(N_m + 1);
    vpz_m.reserve(N_m + 1);

    // reserve memory of N_m for the properties (--> size = 0, capacity = N_m)
    h_m.reserve(N_m);
    ds_m.reserve(N_m);
    fidx_m.reserve(N_m);
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type
    ClosedOrbitFinder<Value_type, Size_type, Stepper>::getInverseBendingRadius(const value_type& angle)
{
    if (angle != 0.0)
        return rotate(angle, h_m);
    else
        return h_m;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type
    ClosedOrbitFinder<Value_type, Size_type, Stepper>::getPathLength(const value_type& angle)
{
    if (angle != 0.0)
        return rotate(angle, ds_m);
    else
        return ds_m;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type
    ClosedOrbitFinder<Value_type, Size_type, Stepper>::getFieldIndex(const value_type& angle)
{
    if (angle != 0.0)
        return rotate(angle, fidx_m);
    return fidx_m;
}

template<typename Value_type, typename Size_type, class Stepper>
std::pair<Value_type,Value_type> ClosedOrbitFinder<Value_type, Size_type, Stepper>::getTunes() {
    // compute radial tune
    value_type nur = computeTune(x_m,px_m[1],nxcross_m);

    // compute vertical tune
    value_type nuz = computeTune(z_m,pz_m[1],nzcross_m);

    return std::make_pair(nur,nuz);
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type
    ClosedOrbitFinder<Value_type, Size_type, Stepper>::getOrbit(value_type angle)
{
    if (angle != 0.0)
        return rotate(angle, r_m);
    else
        return r_m;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type
    ClosedOrbitFinder<Value_type, Size_type, Stepper>::getMomentum(value_type angle)
{
    container_type pr = pr_m;

    if (angle != 0.0)
        pr = rotate(angle, pr);

    // change units from meters to \beta * \gamma
    /* in Gordon paper:
     *
     * p = \gamma * \beta * a
     *
     * where a = c / \omega_{0} with \omega_{0} = 2 * \pi * \nu_{0} = 2 * \pi * \nu_{rf} / h
     *
     * c: speed of light
     * h: harmonic number
     * v_{rf}: nomial rf frequency
     *
     * Units:
     *
     * [a] = m --> [p] = m
     *
     * The momentum in \beta * \gamma is obtained by dividing by "a"
     */
    value_type factor =  1.0 / acon_m(wo_m);
    std::for_each(pr.begin(), pr.end(), [factor](value_type& p) { p *= factor; });

    return pr;
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::value_type
    ClosedOrbitFinder<Value_type, Size_type, Stepper>::getAverageRadius()
{
    return ravg_m;
}

template<typename Value_type, typename Size_type, class Stepper>
typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::value_type
    ClosedOrbitFinder<Value_type, Size_type, Stepper>::getFrequencyError()
{
    return phase_m;
}

// -----------------------------------------------------------------------------------------------------------------------
// PRIVATE MEMBER FUNCTIONS
// -----------------------------------------------------------------------------------------------------------------------



template<typename Value_type, typename Size_type, class Stepper>
bool ClosedOrbitFinder<Value_type, Size_type, Stepper>::findOrbit(value_type accuracy,
                                                                  size_type maxit,
                                                                  value_type ekin,
                                                                  value_type dE,
                                                                  value_type rguess,
                                                                  bool isTuneMode)
{
    /* REMARK TO GORDON
     * q' = 1/b = 1/bcon
     * a' = a = acon
     */

    // resize vectors (--> size = N_m+1, capacity = N_m+1), note: we do N_m+1 integration steps
    r_m.resize(N_m+1);
    pr_m.resize(N_m+1);
    vz_m.resize(N_m+1);
    vpz_m.resize(N_m+1);

    // store acon locally
    value_type acon = acon_m(wo_m);            // [acon] = m
    // amplitude of error; Gordon, formula (18) (a = a')
    value_type error = std::numeric_limits<value_type>::max();

    /*
     * Christian:
     * N = 1440 ---> N = 720 ---> dtheta = 2PI/720 --> nsteps = 721
     *
     * 0, 2, 4, ... ---> jeden zweiten berechnene: 1, 3, 5, ... interpolieren --> 1440 Werte
     *
     * Matthias:
     * N = 1440 --> dtheta = 2PI/1440 --> nsteps = 1441
     *
     * 0, 1, 2, 3, 4, 5, ... --> 1440 Werte
     *
     */

    value_type E;     // starting energy
    value_type E_fin; // final    energy

    if ( isTuneMode ) {
        E     = cycl_m->getFMLowE();
        E_fin = cycl_m->getFMHighE();
    } else {
        E     = ekin;
        E_fin = ekin;
    }

    namespace fs = boost::filesystem;
    fs::path dir = OpalData::getInstance()->getInputBasename();
    dir = dir.parent_path() / "data";
    std::string tunefile = (dir / "tunes.dat").string();

    if ( isTuneMode ) {
        std::ofstream out(tunefile, std::ios::out);

        out << std::left
            << std::setw(15) << "# energy[MeV]"
            << std::setw(15) << "radius_ini[m]"
            << std::setw(15) << "radius_avg[m]"
            << std::setw(15) << "nu_r"
            << std::setw(15) << "nu_z"
            << std::endl;
    }
    // initial guess
    container_type init;
    enum Guess {NONE, FIRST, SECOND};
    Guess guess = NONE;
    value_type rn1, pn1; // normalised r, pr value of previous closed orbit
    value_type rn2, pn2; // normalised r, pr value of second to previous closed orbit

    // iterate until suggested energy (start with minimum energy)
    // increase energy by dE
    for (; E <= E_fin ; E+=dE) {
        error = std::numeric_limits<value_type>::max();

        // energy dependent values
        value_type en     = E / E0_m;        // en = E/E0 = E/(mc^2) (E0 is potential energy)
        value_type gamma  = en + 1.0;
        value_type gamma2 = gamma * gamma;   // = gamma^2
        value_type beta   = std::sqrt(1.0 - 1.0 / gamma2);
        value_type p = acon_m(wo_m) * std::sqrt(en * (2.0 + en));  // momentum [p] = m; Gordon, formula (3)

        if (guess == NONE) {
            // set initial values for radius and radial momentum for lowest energy Emin
            // orbit, [r] = m;  Gordon, formula (20)
            // radial momentum; Gordon, formula (20)
            //      r            pr   z    pz
            init = {beta * acon, 0.0, 0.0, 1.0};
            if (rguess >= 0.0) {
                init[0] = rguess * 0.001;
            }
            guess = FIRST;
        } else if (guess == FIRST) {
            // new initial values based on previous one, formula (21)
            init[0] = (beta*acon) * rn1;
            init[1] = p*pn1;
            guess = SECOND;
        } else if (guess == SECOND) {
            // second extrapolation, formula (21)
            init[0] = (beta*acon) * (rn1 + (rn1-rn2));
            init[1] = p*(pn1 + (pn1-pn2));
        }

        std::fill(  r_m.begin(),   r_m.end(), 0);
        std::fill( pr_m.begin(),  pr_m.end(), 0);
        std::fill( vz_m.begin(),  vz_m.end(), 0);
        std::fill(vpz_m.begin(), vpz_m.end(), 0);

        // (re-)set inital values for r and pr
        r_m[0]   = init[0];
        pr_m[0]  = init[1];
        vz_m[0]  = init[2];
        vpz_m[0] = init[3];

        if ( !this->findOrbitOfEnergy_m(E, init, error, accuracy, maxit) ) {
            *gmsg << "ClosedOrbitFinder didn't converge for energy " + std::to_string(E) + " MeV." << endl;
            guess = NONE;
            continue;
        }
        // store for next initial guess
        rn2 = rn1;
        pn2 = pn1;
        rn1 = r_m[0] / (acon * beta);
        pn1 = pr_m[0] / p;

        if ( isTuneMode ) {

            this->computeOrbitProperties(E);
            value_type reo = this->getOrbit(   cycl_m->getPHIinit())[0];
            value_type peo = this->getMomentum(cycl_m->getPHIinit())[0];
            std::pair<value_type , value_type > tunes = this->getTunes();

            *gmsg << std::left
                  << "* ----------------------------" << endl
                  << "* Closed orbit info (Gordon units):" << endl
                  << "*" << endl
                  << "* kinetic energy:   " << std::setw(12) << E      << " [MeV]" << endl
                  << "* average radius:   " << std::setw(12) << ravg_m << " [m]" << endl
                  << "* initial radius:   " << std::setw(12) << reo    << " [m]" << endl
                  << "* initial momentum: " << std::setw(12) << peo    << " [Beta Gamma]" << endl
                  << "* frequency error:  " << phase_m        << endl
                  << "* horizontal tune:  " << tunes.first    << endl
                  << "* vertical tune:    " << tunes.second   << endl
                  << "* ----------------------------" << endl << endl;

            std::ofstream out(tunefile, std::ios::app);
            out << std::left
                << std::setw(15) << E
                << std::setw(15) << reo
                << std::setw(15) << ravg_m
                << std::setw(15) << tunes.first
                << std::setw(15) << tunes.second << std::endl;
            out.close();
        }
    }

    /* store last entry, since it is needed in computeVerticalOscillations(), because we have to do the same
     * number of integrations steps there.
     */
    /* lastOrbitVal_m    = r_m[N_m];        // needed in computeVerticalOscillations() */
    /* lastMomentumVal_m = pr_m[N_m];       // needed in computeVerticalOscillations() */

    // remove last entry (since we don't have to store [0,2pi], but [0,2pi[)  --> size = N_m, capacity = N_m+1
    r_m.pop_back();
    pr_m.pop_back();

    /* domain_m = true --> only integrated over a single sector
     * --> multiply by cycl_m->getSymmetry() to get correct phase_m
     */
    if (domain_m)
        phase_m *= cycl_m->getSymmetry();

    // returns true if converged, otherwise false
    return error < accuracy;
}

template<typename Value_type, typename Size_type, class Stepper>
bool ClosedOrbitFinder<Value_type, Size_type, Stepper>::findOrbitOfEnergy_m(
    const value_type& E,
    container_type& init,
    value_type& error,
    const value_type& accuracy,
    size_type maxit)
{
    /* *gmsg << "rguess : " << init[0] << endl; */
    /* *gmsg << "prguess: " << init[1] << endl; */

    value_type bint, brint, btint;
    value_type invbcon = 1.0 / bcon_m(E0_m, wo_m);      // [bcon] = MeV*s/(C*m^2) = 10^6 T = 10^7 kG (kilo Gauss)

    value_type xold = 0.0;                              // for counting nxcross
    value_type zold = 0.0;                              // for counting nzcross

    // index for reaching next element of the arrays r and pr (no nicer way found yet)
    size_type idx = 0;
    // observer for storing the current value after each ODE step (e.g. Runge-Kutta step) into the containers of r and pr
    auto store = [&](state_type& y, const value_type t)
    {
        r_m[idx]   = y[0];
        pr_m[idx]  = y[1];
        vz_m[idx]  = y[6];
        vpz_m[idx] = y[7];

        // count number of crossings (excluding starting point --> idx>0)
        nxcross_m += (idx > 0) * (y[4] * xold < 0);
        xold = y[4];

        // number of times z2 changes sign
        nzcross_m += (idx > 0) * (y[10] * zold < 0);
        zold = y[10];

        ++idx;
    };

    // define initial state container for integration: y = {r, pr, x1, px1, x2, px2,
    //                                                      z, pz, z1, pz1, z2, pz2,
    //                                                      phase}
    state_type y(11);

    // difference of last and first value of r (1. element) and pr (2. element)
    container_type err(2);
    // correction term for initial values: r = r + dr, pr = pr + dpr; Gordon, formula (17)
    container_type delta = {0.0, 0.0};
    // if niterations > maxit --> stop iteration
    size_type niterations = 0;

    // energy dependent values
    value_type en     = E / E0_m;                      // en = E/E0 = E/(mc^2) (E0 is potential energy)
    value_type gamma  = en + 1.0;
    value_type p = acon_m(wo_m) * std::sqrt(en * (2.0 + en));     // momentum [p] = m; Gordon, formula (3)
    value_type gamma2    = gamma * gamma;           // = gamma^2
    value_type invgamma4 = 1.0 / (gamma2 * gamma2); // = 1/gamma^4
    value_type p2 = p * p;                              // p^2 = p*p

    // helper constants
    value_type pr2;                                     // squared radial momentum (pr^2 = pr*pr)
    value_type ptheta, invptheta;                       // azimuthal momentum

    // define the six ODEs (using lambda function)
    function_t orbit_integration = [&](const state_type &y,
                                       state_type &dydt,
                                       const double theta)
    {
        pr2 = y[1] * y[1];
        if (p2 < pr2) {
            //*gmsg << theta << " " << p2 << " " << pr2 << endl;
            throw OpalException("ClosedOrbitFinder::findOrbitOfEnergy_m()",
                                "p_{r}^2 > p^{2} (defined in Gordon paper) --> Square root of negative number.");
        }
        // Gordon, formula (5c)
        ptheta    = std::sqrt(p2 - pr2);
        invptheta = 1.0 / ptheta;
        // average field over the number of sectors
        brint=0.0, btint=0.0, bint=0.0;
        for (int i = 0; i<nSectors_m; i++) {
            double angle = theta + i * Physics::two_pi / nSectors_m;
            double tmpbr, tmpbt, tmpb;
            // interpolate values of magnetic field
            cycl_m->apply(y[0], y[6], angle, tmpbr, tmpbt, tmpb);
            brint += tmpbr;
            btint += tmpbt;
            bint  += tmpb;
        }
        brint /= nSectors_m;
        btint /= nSectors_m;
        bint  /= nSectors_m;
        // multiply by 1 / b
        bint  *= invbcon;
        brint *= invbcon;
        btint *= invbcon;

        // Gordon, formula (5a)
        dydt[0] = y[0] * y[1] * invptheta;
        // Gordon, formula (5b) (typo in paper! second equal sign is a minus)
        dydt[1] = ptheta - y[0] * bint;
        // Gordon, formulas (9a) and (9b)
        for (size_type i = 2; i < 5; i += 2) {
            dydt[i]   = (y[1] * y[i] + y[0] * p2 * y[i+1] * invptheta * invptheta) * invptheta;
            dydt[i+1] = - y[1] * y[i+1] * invptheta - (bint + y[0] * brint) * y[i];
        }

        // Gordon, formulas (22a) and (22b)
        for (size_type i = 6; i < 12; i += 2) {
            dydt[i]   = y[0] * y[i+1] * invptheta;
            dydt[i+1] = (y[0] * brint - y[1] * invptheta * btint) * y[i];
        }

        // integrate phase
        dydt[12] = y[0] * invptheta * gamma - 1;

    };

    // integrate until error smaller than user-defined accuracy
    do {
        // (re-)set initial values
        x_m[0]  = 1.0;               // x1;  Gordon, formula (10)
        px_m[0] = 0.0;               // px1; Gordon, formula (10)
        x_m[1]  = 0.0;               // x2;  Gordon, formula (10)
        px_m[1] = 1.0;               // px2; Gordon, formula (10)
        z_m[0]  = 1.0;
        pz_m[0] = 0.0;
        z_m[1]  = 0.0;
        pz_m[1] = 1.0;
        phase_m = 0.0;
        nxcross_m = 0;               // counts the number of crossings of x-axis (excluding first step)
        nzcross_m = 0;
        idx = 0;                     // index for looping over r and pr arrays

        // fill container with initial states
        y = {init[0],init[1],
             x_m[0], px_m[0], x_m[1], px_m[1],
             init[2], init[3],
             z_m[0], pz_m[0], z_m[1], pz_m[1],
             phase_m
        };

        try {
            // integrate from 0 to 2*pi / nSectors (one has to get back to the "origin")
            boost::numeric::odeint::integrate_n_steps(stepper_m, orbit_integration,y,0.0,dtheta_m,N_m,store);
        } catch(OpalException & ex) {
            *gmsg << ex.where() << " " << ex.what() << endl;
            break;
        }

        // write new state
        x_m[0]  = y[2];
        px_m[0] = y[3];
        x_m[1]  = y[4];
        px_m[1] = y[5];

        z_m[0]  = y[8];
        pz_m[0] = y[9];
        z_m[1]  = y[10];
        pz_m[1] = y[11];
        phase_m = y[12] * Physics::u_two_pi; // / (2.0 * Physics::pi);

        // compute error (compare values of orbit and momentum for theta = 0 and theta = 2*pi)
        // (Note: size = N_m+1 --> last entry is N_m)
        err[0] =  r_m[N_m] -  r_m[0];    // Gordon, formula (14)
        err[1] = pr_m[N_m] - pr_m[0];    // Gordon, formula (14)

        // correct initial values of r and pr
        value_type invdenom = 1.0 / (x_m[0] + px_m[1] - 2.0);
        delta[0] = ((px_m[1] - 1.0) * err[0] -  x_m[1] * err[1]) * invdenom; // dr;  Gordon, formula (16a)
        delta[1] = (( x_m[0] - 1.0) * err[1] - px_m[0] * err[0]) * invdenom; // dpr; Gordon, formula (16b)

        // improved initial values; Gordon, formula (17) (here it's used for higher energies)
        init[0] += delta[0];
        init[1] += delta[1];

        // compute amplitude of the error (Gordon, formula (18)
        error = std::sqrt(delta[0] * delta[0] + delta[1] * delta[1] * invgamma4) / r_m[0];

        // *gmsg << "iteration " << niterations << " error: " << error << endl;

    } while ((error > accuracy) && (niterations++ < maxit));

    if (error > accuracy)
        *gmsg << "findOrbit not converged after " << niterations << " iterations with error: " << error << ". Needed accuracy " << accuracy << endl;

    return (error < accuracy);
}

template<typename Value_type, typename Size_type, class Stepper>
Value_type ClosedOrbitFinder<Value_type, Size_type, Stepper>::computeTune(const std::array<value_type,2>& y,
                                                                          value_type py2, size_type ncross)
{
    // Y = [y1, y2; py1, py2]

    // cos(mu)
    value_type cos = 0.5 * (y[0] + py2);

    value_type mu;

    // sign of sin(mu) has to be equal to y2
    bool negative = std::signbit(y[1]);

    bool uneven = (ncross % 2);

    value_type abscos = std::fabs(cos);
    value_type muPrime;
    if (abscos > 1.0) {
        // store the number of crossings
        if (uneven)
            ncross = ncross + 1;

        // Gordon, formula (36b)
        muPrime = -std::acosh(abscos);    // mu'

    } else {
        muPrime = (uneven) ? std::acos(-cos) : std::acos(cos);    // mu'
        /* It has to be fulfilled: 0<= mu' <= pi
         * But since |cos(mu)| <= 1, we have
         * -1 <= cos(mu) <= 1 --> 0 <= mu <= pi (using above programmed line), such
         * that condition is already fulfilled.
         */
    }

    // Gordon, formula (36)
    mu = ncross * Physics::pi + muPrime;

    if (abscos < 1.0) {
        // if sign(y[1]) > 0 && sign(sin(mu)) < 0
        if (!negative && std::signbit(std::sin(mu))) {
            mu = ncross * Physics::pi - muPrime;
        } else if (negative && !std::signbit(std::sin(mu))) {
            mu = ncross * Physics::pi - muPrime + Physics::two_pi;
        }
    }

    // nu = mu/theta, where theta = integration domain

    /* domain_m = true --> only integrated over a single sector --> multiply by cycl_m->getSymmetry() to
     * get correct tune.
     */
    if (domain_m)
        mu *= cycl_m->getSymmetry();

    return mu * Physics::u_two_pi;
}

template<typename Value_type, typename Size_type, class Stepper>
void ClosedOrbitFinder<Value_type, Size_type, Stepper>::computeOrbitProperties(const value_type& E) {
    /*
     * The formulas for h, fidx and ds are from the paper:
     * "Tranverse-Longitudinal Coupling by Space Charge in Cyclotrons"
     * written by Dr. Christian Baumgarten (2012, PSI)
     * p. 6
     */

    // READ IN MAGNETIC FIELD: ONLY FOR STAND-ALONE PROGRAM
    value_type bint, brint, btint; // B, dB/dr, dB/dtheta

    value_type invbcon = 1.0 / bcon_m(E0_m, wo_m);
    value_type en = E / E0_m;                                  // en = E/E0 = E/(mc^2) (E0 is potential energy)
    value_type p = acon_m(wo_m) * std::sqrt(en * (2.0 + en));  // momentum [p] = m; Gordon, formula (3)
    value_type p2 = p * p;
    value_type theta = 0.0;                                    // angle for interpolating
    value_type ptheta;

    // resize of container (--> size = N_m, capacity = N_m)
    h_m.resize(N_m);
    fidx_m.resize(N_m);
    ds_m.resize(N_m);

    for (size_type i = 0; i < N_m; ++i) {
        // interpolate magnetic field
        cycl_m->apply(r_m[i], vz_m[i], theta, brint, btint, bint);
        bint  *= invbcon;
        brint *= invbcon;
        btint *= invbcon;

        // inverse bending radius
        h_m[i] = bint / p;

        // local field index
        if (p < pr_m[i])
            throw OpalException("ClosedOrbitFinder::computeOrbitProperties()",
                                "p_{r}^2 > p^{2} " + std::to_string(p) + " " + std::to_string(pr_m[i]) + " (defined in Gordon paper) --> Square root of negative number.");

        ptheta = std::sqrt(p2 - pr_m[i] * pr_m[i]);

        fidx_m[i] = (brint * ptheta - btint * pr_m[i] / r_m[i]) / p2; //(bint*bint);

        // path length element
        ds_m[i] = std::hypot(r_m[i] * pr_m[i] / ptheta,r_m[i]) * dtheta_m; // C++11 function

        // increase angle
        theta += dtheta_m;
    }

    // compute average radius
    ravg_m = std::accumulate(r_m.begin(),r_m.end(),0.0) / value_type(r_m.size());
}

template<typename Value_type, typename Size_type, class Stepper>
inline typename ClosedOrbitFinder<Value_type, Size_type, Stepper>::container_type
ClosedOrbitFinder<Value_type, Size_type, Stepper>::rotate(value_type angle, container_type &orbitProperty) {

    container_type orbitPropertyCopy = orbitProperty;

    // compute the number of steps per degree
    value_type deg_step = N_m / 360.0;

    // compute starting point
    size_type start = deg_step * angle;

    // copy end to start
    std::copy(orbitProperty.begin() + start, orbitProperty.end(), orbitPropertyCopy.begin());

    // copy start to end
    std::copy_n(orbitProperty.begin(), start, orbitPropertyCopy.end() - start);

    return orbitPropertyCopy;

}

#endif
