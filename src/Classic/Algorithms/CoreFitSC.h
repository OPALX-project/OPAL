//
// Namespace CoreFitSC
//   Core selection and far-field model for a space-charge mesh fitted to the bunch core
//   (OPAL-cycl, FFT solver).
//
//   The mesh is fitted to the core: the live particles within nSigma clipped rms widths of
//   the clipped mean in every axis of the local space-charge frame. The clipped mean and rms
//   come from an iterated clip at clip rms widths that starts from the plain rms, so a few
//   particles far out in the halo neither stretch the mesh nor inflate its cell size. The
//   other particles ("far") are not deposited. Within the range of the mesh's cell centres
//   they are gathered from the mesh, as the core is. Outside it they get the rest-frame
//   monopole and quadrupole field of the core charge, or, where that expansion is not valid,
//   a direct sum over the core particles. Every far particle also gets a direct sum over the
//   other live far particles. The lab-frame transform is the one of the mesh field.
//
//   Everything here is a pure function of its arguments. Global sums and minima go through
//   an injected Reducer, so the same code runs serially in the unit tests and on all ranks
//   in OPAL. Nothing is kept between calls: every selection starts from the plain rms, and
//   the result depends on the particle set only (and on the rank count only at round-off).
//   Positions are in the local space-charge frame in m, charges in C, and fields are
//   rest-frame fields without the coupling constant, sum q d / |d|^3 in C/m^2.
//
// Created by: Daniel Winklehner, MIT (2026)
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#ifndef CLASSIC_COREFITSC_H
#define CLASSIC_COREFITSC_H

#include "Algorithms/Vektor.h"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace CoreFitSC {

    /// Fixed parameters of the selection and of the far-field model.
    /// At most this many clipped passes after the plain rms.
    constexpr unsigned int maxClipPasses = 8;
    /// The clipped passes stop when |dsigma_d| < clipTolerance * sigma_d in every axis.
    constexpr double clipTolerance = 1e-3;
    /// The halo proper: live particles beyond haloSigma clipped rms widths in any axis.
    constexpr double haloSigma = 10.0;
    /// Smallest half-width of the clip, core and halo tests [m]. In an axis in which the bunch
    /// is flat (rms 0 up to round-off), the tests then keep the particles in its plane,
    /// whatever the round-off of the mean, and exclude the others. boundp() treats extents
    /// below the same 1e-10 m as flat; the rms of a real bunch is many orders larger.
    constexpr double flatWidth = 1e-10;
    /// FAR particles whose quadrupole term exceeds this fraction of the monopole term get
    /// the exact core field (needsExactCoreField()).
    constexpr double indicatorThreshold = 0.1;

    /// Mode decision: the solve falls back to the full box (FULL) when the core holds
    /// fewer than max(minCoreParticles, minCoreFraction * live) particles or the far charge
    /// fraction exceeds its limit. It returns to the core box (CORE) only when the far charge
    /// fraction is below returnFarFactor times the limit and the core holds more than
    /// returnCoreFraction of the live particles (and at least minCoreParticles).
    constexpr double minCoreParticles = 1000;
    constexpr double minCoreFraction = 0.5;
    constexpr double returnCoreFraction = 0.6;
    constexpr double returnFarFactor = 0.5;

    /// Global reductions in place over all ranks. Every rank must make the same calls with
    /// the same counts and must get back the same values, as MPI_Allreduce returns them.
    /// The selection takes every decision from reduced values only, so the ranks stay in
    /// step.
    class Reducer {
    public:
        virtual ~Reducer() = default;
        /// Sums over all ranks
        virtual void sum(double* values, unsigned int count) const = 0;
        /// Minima over all ranks (maxima are the minima of the negated values)
        virtual void min(double* values, unsigned int count) const = 0;
    };

    /// All particles on one process: the reductions do nothing.
    class SerialReducer: public Reducer {
    public:
        void sum(double*, unsigned int) const override { }
        void min(double*, unsigned int) const override { }
    };

    /// All IPPL nodes, through allreduce().
    class IpplReducer: public Reducer {
    public:
        void sum(double* values, unsigned int count) const override;
        void min(double* values, unsigned int count) const override;
    };

    /// Result of selectCore(), the same on every rank.
    struct CoreSelection {
        Vector_t mean;             ///< clipped mean [m]
        Vector_t sigma;            ///< clipped rms [m]
        Vector_t boundsMin;        ///< minimum over the core particles [m]
        Vector_t boundsMax;        ///< maximum over the core particles [m]; an empty box
                                   ///< (boundsMin > boundsMax) if there is no core particle
        std::size_t numLive = 0;   ///< particles with Bin >= 0
        std::size_t numCore = 0;   ///< live particles within nSigma * sigma in every axis
        std::size_t numFar = 0;    ///< the other live particles
        std::size_t numHalo = 0;   ///< live particles beyond haloSigma * sigma in any axis
        double chargeLive = 0.0;   ///< charge of the live particles [C]
        double chargeFar = 0.0;    ///< charge of the far live particles [C]
        double farMaxDistance = 0.0;  ///< largest max_d |r_d - mean_d| / sigma_d of a far
                                      ///< particle, over the axes with sigma_d > flatWidth
        unsigned int passes = 0;   ///< clipped passes done, 0 to maxClipPasses
        bool converged = false;    ///< false if the passes stopped at maxClipPasses

        /// Far charge over live charge, 0 without live charge
        double getFarChargeFraction() const;
    };

    /// Core selection over the localNum particles of this rank (local space-charge frame).
    /// Pass 0 is the plain mean and rms over the live particles (Bin >= 0). Each of the
    /// passes 1 to maxClipPasses takes the mean and rms over the live particles within
    /// clip * sigma of the mean of the previous pass in every axis. The passes stop when the
    /// count does not change or when every sigma_d changed by less than clipTolerance (an
    /// axis whose new sigma_d is at most flatWidth has settled); a pass that selects no
    /// particle is dropped and ends them. The core is the live particles within
    /// nSigma * sigma in every axis, and its bounds are the exact minima and maxima of their
    /// positions. The halo counts every live particle beyond haloSigma * sigma in any axis,
    /// core particles too if nSigma > haloSigma. No half-width of these tests is below
    /// flatWidth. Particles flagged lost count nowhere. The pointers may be null if localNum
    /// is 0. Collectives: 2 + passes sums (3 + passes if a pass is dropped) and one minimum,
    /// of at most 8 doubles.
    CoreSelection selectCore(const Vector_t* R, const double* Q, const int* bin,
                             std::size_t localNum, double nSigma, double clip,
                             const Reducer& reducer);

    /// True if r lies in [lower, upper] in every axis
    inline bool isInside(const Vector_t& r, const Vector_t& lower, const Vector_t& upper) {
        return (r[0] >= lower[0] && r[0] <= upper[0] &&
                r[1] >= lower[1] && r[1] <= upper[1] &&
                r[2] >= lower[2] && r[2] <= upper[2]);
    }

    /// How a particle takes part in the solve on the core mesh. The far table holds the
    /// MESH and FAR particles, with charge 0 for the lost ones.
    enum class ParticleClass: short {
        CORE,     ///< live, inside the core bounds: deposited, field gathered from the mesh
        INSIDE,   ///< lost (Bin < 0), inside the core bounds: gathered, not deposited
        MESH,     ///< outside the core bounds, inside the range [meshMin, meshMax] of the
                  ///< mesh's cell centres: gathered, plus the far-far field
        FAR       ///< outside that range: multipole (or exact core) field, plus far-far
    };

    /// Class of a particle, from the core bounds of the selection and the cell-centre range
    /// of the core mesh (the core bounds enlarged by boundp()). A live particle is CORE
    /// exactly when selectCore() put it in the core: a far particle fails the nSigma test in
    /// some axis, and the bounds in that axis are extremes of particles that pass it.
    ParticleClass classify(int bin, const Vector_t& r,
                           const Vector_t& boundsMin, const Vector_t& boundsMax,
                           const Vector_t& meshMin, const Vector_t& meshMax);

    /// Rest-frame position about the centre r0 (the core box centre),
    /// s = (x - x0, gamma (y - y0), z - z0): the frame of the mesh solve, whose longitudinal
    /// mesh spacing is hr[1] * gamma (PartBunch::computeSelfFields_cycl).
    inline Vector_t toRestFrame(const Vector_t& r, const Vector_t& r0, double gamma) {
        return Vector_t({r[0] - r0[0], gamma * (r[1] - r0[1]), r[2] - r0[2]});
    }

    /// Lab-frame field of a rest-frame field, (gamma E'x, E'y, gamma E'z): the factors the
    /// mesh field gets from eg_m *= (gamma, 1/gamma, gamma) after the gradient over the
    /// unscaled mesh spacing.
    inline Vector_t toLabFrame(const Vector_t& E, double gamma) {
        return Vector_t({gamma * E[0], E[1], gamma * E[2]});
    }

    /// Rest-frame moments of the core charge.
    struct CoreMoments {
        double charge = 0.0;                  ///< Q_c [C]
        Vector_t centroid;                    ///< c', rest frame, about the centre [m]
        std::array<Vector_t, 3> quadrupole;   ///< Q_ij = 3 M_ij - tr(M) delta_ij about c',
                                              ///< M_ij = sum q s_i s_j - Q_c c'_i c'_j [C m^2]
    };

    /// Local sums of q, q s and q s_i s_j (10 numbers) over the particles added, with s the
    /// rest-frame position about a fixed centre (the core box centre, for conditioning).
    class MomentSums {
    public:
        MomentSums(const Vector_t& centre, double gamma);

        /// Adds a particle of charge q at r (local frame, not rest frame)
        void add(double q, const Vector_t& r);

        /// Moments of all ranks' particles: one sum of 10 doubles
        CoreMoments reduce(const Reducer& reducer) const;

    private:
        Vector_t centre_m;
        double gamma_m;
        std::array<double, 10> sums_m;
    };

    /// Monopole plus quadrupole field of the core at the rest-frame positions of the far
    /// table entries [begin, end), without the coupling constant:
    ///   E' = Q_c d / |d|^3 - Q d / |d|^5 + (5/2) (d Q d) d / |d|^7,  d = s - c',
    /// minus the gradient of the potential Q_c / |d| + (1/2) d Q d / |d|^5 (the dipole about
    /// c' is zero). indicator[i] is |quadrupole term| / |monopole term| (0 without monopole).
    /// With withQuadrupole false, field[i] is the monopole term alone. field and indicator
    /// must have at least end entries; entries outside [begin, end) are not touched. Only the
    /// FAR entries use them: a FAR entry with indicator > indicatorThreshold gets the exact
    /// core field instead (needsExactCoreField()), and the MESH entries are gathered from the
    /// mesh.
    void computeMultipoleField(const CoreMoments& moments, const std::vector<Vector_t>& position,
                               std::size_t begin, std::size_t end, bool withQuadrupole,
                               std::vector<Vector_t>& field, std::vector<double>& indicator);

    /// True if a far table entry gets the exact core field in place of the multipole field:
    /// a FAR entry whose indicator exceeds indicatorThreshold. A MESH entry never does: it is
    /// gathered from the mesh, whatever its indicator.
    inline bool needsExactCoreField(ParticleClass c, double indicator) {
        return c == ParticleClass::FAR && indicator > indicatorThreshold;
    }

    /// Adds to field[i], for the targets i in [begin, end), the rest-frame field of the point
    /// charges at source, without the coupling constant, softened by eps:
    ///   sum_j charge_j (t_i - s_j) / (|t_i - s_j|^2 + eps^2)^(3/2).
    /// Pairs at exactly the same position are skipped. With the far table as targets and its
    /// live entries as sources this is the far-far field (eps = the smallest rest-frame core
    /// cell); the lost entries, of charge 0, would only add cost as sources. With the local
    /// core particles as sources and eps = 0 it is this rank's part of the exact core field.
    void addDirectField(const std::vector<Vector_t>& target, std::size_t begin, std::size_t end,
                        const std::vector<Vector_t>& source, const std::vector<double>& charge,
                        double eps, std::vector<Vector_t>& field);

    /// Entries [rank * n / numRanks, (rank + 1) * n / numRanks) of a table of n entries: the
    /// far-field work of one rank, independent of which rank owns the particles.
    std::pair<std::size_t, std::size_t> getTableRange(std::size_t n, int rank, int numRanks);

    /// CORE solves on the core-fitted mesh, FULL on the box of all particles (fallback).
    /// The values are those of the stat column scMode.
    enum class Mode: short {
        FULL = 0,
        CORE = 1
    };

    /// Mode of the solve and the steps solved in it since the last switch, the switching
    /// step included; dwell saturates at max(minDwell, 1). dwell 0 means that the run has
    /// not solved yet and counts as the minimum dwell reached: a run starts from
    /// ModeState{} = {CORE, 0}, and its first solve may switch.
    struct ModeState {
        Mode mode = Mode::CORE;
        unsigned int dwell = 0;
    };

    /// Mode for the next solve, from globally reduced inputs and the previous state, so every
    /// rank takes the same decision. CORE switches to FULL at once when numCore <
    /// max(minCoreParticles, minCoreFraction * numLive): the core has broken up, or it is
    /// empty and has no box. The other switches wait until dwell reaches minDwell
    /// (REPARTFREQ): CORE to FULL when farChargeFraction > maxFarFraction, and FULL back to
    /// CORE when farChargeFraction < returnFarFactor * maxFarFraction, numCore >
    /// returnCoreFraction * numLive and numCore >= minCoreParticles. The mode switched if it
    /// differs from previous.mode.
    ModeState decideMode(const ModeState& previous, std::size_t numCore, std::size_t numLive,
                         double farChargeFraction, double maxFarFraction, unsigned int minDwell);
}

#endif
