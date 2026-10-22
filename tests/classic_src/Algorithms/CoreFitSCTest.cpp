//
// Unit tests for CoreFitSC, the core selection and far-field model of the core-fitted
// space-charge mesh.
//
//   The expected values of the synthetic bunches were printed by an independent numpy
//   implementation of the same rule and model, which is not part of OPAL. Both generate
//   the bunches with the same splitmix64 / Box-Muller generator, so the particles are
//   bit-identical, and the results must agree to round-off. As OPAL does, the numpy
//   implementation flags only the FAR entries of the far table for the exact core field;
//   the MESH entries are gathered from the mesh.
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
#include "gtest/gtest.h"

#include "Algorithms/CoreFitSC.h"
#include "Physics/Physics.h"
#include "Utility/IpplInfo.h"

#include <algorithm>
#include <atomic>
#include <barrier>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <thread>
#include <vector>

using namespace CoreFitSC;

namespace {
    /// Output function of splitmix64
    std::uint64_t mix64(std::uint64_t z) {
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }

    /// Hash of a list of indices, as the numpy implementation computes it
    std::uint64_t hashIndices(const std::vector<std::size_t>& values) {
        std::uint64_t h = 0;
        for (std::size_t v: values) {
            h = mix64(h ^ (v + 0x9E3779B97F4A7C15ULL));
        }
        return h;
    }

    /// splitmix64 with Box-Muller normals, the generator of the numpy implementation
    class SplitMix64 {
    public:
        explicit SplitMix64(std::uint64_t seed): state_m(seed) { }

        std::uint64_t next() {
            state_m += 0x9E3779B97F4A7C15ULL;
            return mix64(state_m);
        }

        double uniform() {
            return (next() >> 11) * (1.0 / 9007199254740992.0);
        }

        double normal() {
            if (hasSpare_m) {
                hasSpare_m = false;
                return spare_m;
            }
            const double u1 = uniform();
            const double u2 = uniform();
            const double rad = std::sqrt(-2.0 * std::log(1.0 - u1));
            const double phi = 2.0 * Physics::pi * u2;
            spare_m = rad * std::sin(phi);
            hasSpare_m = true;
            return rad * std::cos(phi);
        }

    private:
        std::uint64_t state_m;
        double spare_m = 0.0;
        bool hasSpare_m = false;
    };

    struct Bunch {
        std::vector<Vector_t> R;
        std::vector<double> Q;
        std::vector<int> bin;

        void push(const Vector_t& r, double q, int b) {
            R.push_back(r);
            Q.push_back(q);
            bin.push_back(b);
        }
    };

    const Vector_t syntheticSigma({1.0e-3, 3.0e-3, 2.5e-3});
    const double syntheticCharge = 1.0e-15;
    const double syntheticGamma = 1.0641;
    const std::size_t syntheticCore = 200000;

    /// The Gaussian (sigma 1, 3, 2.5 mm) plus 40 halo particles at 15-35 sigma, or the
    /// scraped variant (x truncated at +1.5 sigma, plus 1 % tails 4x wider); then 10
    /// particles flagged lost at 20-50 sigma.
    Bunch makeSynthetic(bool scraped) {
        SplitMix64 g(20260924 + (scraped ? 1 : 0));
        const Vector_t& sig = syntheticSigma;
        Bunch b;
        while (b.R.size() < syntheticCore) {
            const double x = g.normal() * sig[0];
            const double y = g.normal() * sig[1];
            const double z = g.normal() * sig[2];
            if (scraped && x > 1.5 * sig[0]) {
                continue;
            }
            b.push(Vector_t({x, y, z}), syntheticCharge, 0);
        }
        // radius log-uniform in [lo, hi); logRatio = log(hi / lo) is a literal, as in the
        // numpy implementation, so that no compile-time folding can change a bit
        auto shell = [&g, &sig](double lo, double logRatio) {
            const double u0 = g.normal();
            const double u1 = g.normal();
            const double u2 = g.normal();
            const double m = std::max({std::abs(u0), std::abs(u1), std::abs(u2)});
            const double rho = lo * std::exp(g.uniform() * logRatio);
            return Vector_t({rho * u0 / m * sig[0], rho * u1 / m * sig[1], rho * u2 / m * sig[2]});
        };
        if (scraped) {
            for (std::size_t i = 0; i < syntheticCore / 100; ++i) {
                const double x = 4.0 * g.normal() * sig[0];
                const double y = 4.0 * g.normal() * sig[1];
                const double z = 4.0 * g.normal() * sig[2];
                b.push(Vector_t({x, y, z}), syntheticCharge, 0);
            }
        } else {
            for (int i = 0; i < 40; ++i) {
                b.push(shell(15.0, 0.8472978603872037), syntheticCharge, 0);
            }
        }
        for (int i = 0; i < 10; ++i) {
            b.push(shell(20.0, 0.9162907318741551), syntheticCharge, -1);
        }
        return b;
    }

    /// 100000 Gaussian particles, where for the given fraction x is a Cauchy variable (a
    /// ratio of normals) times sigma_x
    Bunch makeHeavyTailed(double fraction) {
        SplitMix64 g(20260926);
        const Vector_t& sig = syntheticSigma;
        Bunch b;
        for (int i = 0; i < 100000; ++i) {
            double x = g.normal() * sig[0];
            const double y = g.normal() * sig[1];
            const double z = g.normal() * sig[2];
            if (g.uniform() < fraction) {
                const double numerator = g.normal();
                x = (numerator / g.normal()) * sig[0];
            }
            b.push(Vector_t({x, y, z}), syntheticCharge, 0);
        }
        return b;
    }

    /// The Gaussian bunch flat in z, z = z0 for every particle, or with coreOnly for the
    /// core particles only, while the others move to z + z0
    Bunch makeFlat(const Bunch& gaussian, double z0, bool coreOnly) {
        Bunch b = gaussian;
        for (std::size_t i = 0; i < b.R.size(); ++i) {
            if (coreOnly && i >= syntheticCore) {
                b.R[i][2] += z0;
            } else {
                b.R[i][2] = z0;
            }
        }
        return b;
    }

    CoreSelection select(const Bunch& b, const Reducer& reducer = SerialReducer(),
                         double nSigma = 6.0) {
        return selectCore(b.R.data(), b.Q.data(), b.bin.data(), b.R.size(), nSigma, 4.0,
                          reducer);
    }

    /// Moments of the live particles inside the core bounds, about the box centre
    CoreMoments getCoreMoments(const Bunch& b, const CoreSelection& s, const Reducer& reducer) {
        MomentSums sums(0.5 * (s.boundsMin + s.boundsMax), syntheticGamma);
        for (std::size_t i = 0; i < b.R.size(); ++i) {
            if (b.bin[i] >= 0 && isInside(b.R[i], s.boundsMin, s.boundsMax)) {
                sums.add(b.Q[i], b.R[i]);
            }
        }
        return sums.reduce(reducer);
    }

    /// The far-field model on a whole bunch, serially, as in the numpy implementation: mesh
    /// cell-centre range as boundp() enlarges the core bounds (BBOXINCR 2 %, 36 cells),
    /// moments of the core about the box centre, far table in particle order with its live
    /// entries as far-far sources, and the FAR entries flagged by the indicator.
    struct FarModel {
        CoreSelection selection;
        Vector_t meshMin, meshMax, hr, centre;
        CoreMoments moments;
        std::vector<ParticleClass> classes;
        std::vector<std::size_t> table;
        std::size_t numMesh = 0;
        std::vector<Vector_t> position;
        std::vector<double> charge;
        std::vector<Vector_t> farSource;
        std::vector<double> farCharge;
        std::vector<Vector_t> coreSource;
        std::vector<double> coreCharge;
        std::vector<Vector_t> multipole;
        std::vector<double> indicator;
        std::vector<Vector_t> farFar;
        std::vector<std::size_t> flagged;   // table rows
        double eps = 0.0;
    };

    FarModel computeFarModel(const Bunch& b, double gamma) {
        FarModel m;
        m.selection = select(b);
        const double dh = 0.02;
        const int nr = 36;
        for (unsigned int d = 0; d < 3; ++d) {
            const double length = std::abs(m.selection.boundsMax[d] - m.selection.boundsMin[d]);
            m.meshMin[d] = m.selection.boundsMin[d] - dh * length;
            m.meshMax[d] = m.selection.boundsMax[d] + dh * length;
            m.hr[d] = (m.meshMax[d] - m.meshMin[d]) / (nr - 1);
        }
        m.centre = 0.5 * (m.selection.boundsMin + m.selection.boundsMax);
        m.eps = std::min({m.hr[0], gamma * m.hr[1], m.hr[2]});

        MomentSums sums(m.centre, gamma);
        for (std::size_t i = 0; i < b.R.size(); ++i) {
            const ParticleClass c = classify(b.bin[i], b.R[i],
                                             m.selection.boundsMin, m.selection.boundsMax,
                                             m.meshMin, m.meshMax);
            m.classes.push_back(c);
            const Vector_t s = toRestFrame(b.R[i], m.centre, gamma);
            if (c == ParticleClass::CORE) {
                sums.add(b.Q[i], b.R[i]);
                m.coreSource.push_back(s);
                m.coreCharge.push_back(b.Q[i]);
            } else if (c == ParticleClass::MESH || c == ParticleClass::FAR) {
                m.table.push_back(i);
                m.position.push_back(s);
                m.charge.push_back(b.bin[i] >= 0 ? b.Q[i] : 0.0);
                if (b.bin[i] >= 0) {
                    m.farSource.push_back(s);
                    m.farCharge.push_back(b.Q[i]);
                }
                if (c == ParticleClass::MESH) {
                    ++ m.numMesh;
                }
            }
        }
        m.moments = sums.reduce(SerialReducer());

        const std::size_t n = m.table.size();
        m.multipole.assign(n, Vector_t(0.0));
        m.indicator.assign(n, 0.0);
        m.farFar.assign(n, Vector_t(0.0));
        computeMultipoleField(m.moments, m.position, 0, n, true, m.multipole, m.indicator);
        addDirectField(m.position, 0, n, m.farSource, m.farCharge, m.eps, m.farFar);
        for (std::size_t k = 0; k < n; ++k) {
            if (needsExactCoreField(m.classes[m.table[k]], m.indicator[k])) {
                m.flagged.push_back(k);
            }
        }
        return m;
    }

    /// Selection values printed by the numpy implementation
    struct ExpectedSelection {
        unsigned int passes;
        bool converged;
        std::size_t numLive, numCore, numFar, numHalo;
        double chargeLive, chargeFar, farMaxDistance;
        double mean[3], sigma[3], boundsMin[3], boundsMax[3];
    };

    /// Selection and model values printed by the numpy implementation
    struct Expected {
        ExpectedSelection selection;
        std::size_t tableSize, numMesh, numFlagged;
        std::uint64_t tableHash;    // hashIndices(4 * particle index + class), table order
        std::uint64_t flaggedHash;  // hashIndices(flagged table rows)
        double coreCharge;
        double centroid[3];
        double quadrupole[6];   // xx, xy, xz, yy, yz, zz
        double eps;
        std::size_t rows[4];    // table rows: first, second, middle, last
        double multipole[4][3];
        double indicator[4];
        double farFar[4][3];
        double coreField[4][3]; // exact core field (eps = 0) at the rows
        double sums[4];         // sum |E_mp|, sum |E_ff|, sum eta, max eta over the table
    };

    // table 50 = 0 MESH + 50 FAR entries; 0 FAR entries flagged
    const Expected gaussianExpected = {
        {2, true, 200040, 200000, 40, 40, 2.0003999999999996e-10, 4.0000000000000006e-14,
         34.89918019251822, {3.622343297332285e-07, 1.534927388850182e-06, 5.617058889432891e-06},
         {0.0009959923186301967, 0.003000495608449616, 0.0024952342662769283},
         {-0.004730037975479777, -0.013018094678226786, -0.01138603316183725},
         {0.004442721171545147, 0.013430044993129886, 0.011973815356987829}},
        50, 0, 0, 0x98de47fe40058803ULL, 0x0000000000000000ULL,
        1.9999999999999996e-10,
        {0.00014393929294318386, -0.0002174979073846089, -0.000288303482864086},
        {-2.890441600205627e-15, -2.2864316108728414e-18, 5.081183124686262e-18,
         2.6385550851631385e-15, 3.3109013006512145e-18, 2.518865150424885e-16},
        0.0002725619860830263,
        {0, 1, 25, 49},
        {{2.5699536046509824e-08, -7.73207535250013e-08, -1.6342013065766556e-08},
         {1.198134310202825e-08, 3.265204136237126e-08, -4.649345374931572e-09},
         {-3.3526039279383695e-09, -6.295938637987614e-08, -2.421316458351747e-08},
         {1.0840939359798923e-08, 2.7587216293127342e-08, 1.9820224625830157e-08}},
        {0.007147827812336598, 0.003003017490984653, 0.006013702066418084, 0.0024224855174186455},
        {{5.493116632237291e-12, 5.972210064125615e-12, -5.2520738143312826e-12},
         {3.8478657438571934e-12, 1.2589983856917994e-12, 1.5253812570167292e-12},
         {4.346097475788938e-12, -1.2311378365940645e-12, -8.717081082907719e-13},
         {2.7630452448754765e-12, 3.1338388268793964e-12, 2.589010362795313e-12}},
        {{2.570732547564313e-08, -7.732138746544804e-08, -1.6343191538024863e-08},
         {1.1981917188110592e-08, 3.2651990960254204e-08, -4.649341457392467e-09},
         {-3.3535381028272597e-09, -6.296288187648848e-08, -2.4215882137002324e-08},
         {1.0841267741934942e-08, 2.7587035125138585e-08, 1.982019350506708e-08}},
        {2.8186423360157568e-06, 2.878470454072526e-10, 0.20031000868695767, 0.03146190802108764}
    };

    // table 769 = 62 MESH + 707 FAR entries; 90 FAR entries flagged, while the 9 MESH
    // entries above the threshold are gathered from the mesh
    const Expected scrapedExpected = {
        {3, true, 202000, 201241, 759, 95, 2.0199999999999997e-10, 7.590000000000002e-13,
         14.560080083496405,
         {-0.00013782758998007722, -6.16470676303518e-06, 1.4220842864287944e-06},
         {0.0008847595290983996, 0.0030115993205686808, 0.0025127535846909146},
         {-0.005444287642411073, -0.018049872500206163, -0.015004148476957541},
         {0.005157099722333585, 0.01803308600152765, 0.0150250846180197}},
        769, 62, 90, 0xef74cd957caadc48ULL, 0xe6b725ad5d909eafULL,
        2.0124099999999998e-10,
        {5.627224291397971e-06, 4.203558947571086e-06, -1.208864107780088e-05},
        {-3.1313219804604118e-15, 3.8848883608914863e-19, -5.5189553136878204e-18,
         2.804113089152963e-15, -6.23905758205239e-18, 3.2720889130744884e-16},
        0.0003150126531238413,
        {0, 1, 384, 768},
        {{3.029403446679831e-08, 1.7497452022624514e-07, 5.673535388574859e-07},
         {9.679071773818593e-08, -5.298876686817554e-08, 3.277470982479987e-07},
         {1.5204936938216584e-07, -2.119251593877244e-07, 3.015901162467223e-07},
         {-1.0659834652019833e-09, -2.993362226915319e-08, 2.8040599709835715e-09}},
        {0.015939618317048248, 0.008744216591131034, 0.019934901588826548, 0.003092488002536406},
        {{2.5549779532506316e-09, 7.239103517535172e-11, 1.057102979226418e-09},
         {6.754883782160742e-10, -2.402388803735614e-10, 1.0417906757036184e-09},
         {8.601610088296822e-10, -4.305341484202658e-10, 6.415764105439749e-10},
         {-5.816752618294435e-12, -1.2285192922370683e-10, 1.2044639086135127e-11}},
        {{3.039259895519485e-08, 1.752203070763518e-07, 5.678117411112008e-07},
         {9.691924434796106e-08, -5.2994033747174154e-08, 3.277743928775997e-07},
         {1.5236171890234047e-07, -2.1175990702924775e-07, 3.0127015949944524e-07},
         {-1.0661220211887495e-09, -2.993397703982991e-08, 2.8040855181848964e-09}},
        {0.0005518596749524464, 8.383812538180716e-07, 42.34560885434678, 0.5842537460402216}
    };

    // heavy tails, fraction 0.25: passes 8 converged 1, core 97650 far 2350 halo 1416
    const ExpectedSelection heavy25Expected =
        {8, true, 100000, 97650, 2350, 1416, 1.0000000000000003e-10, 2.3500000000000007e-12,
         10996.957287674104,
         {-4.716932628771225e-06, 9.830378358713737e-06, -1.6699398414101466e-06},
         {0.0011405297184202334, 0.002992356800897916, 0.0025070459742445694},
         {-0.006843561554146823, -0.012445697203554977, -0.00977668221053564},
         {0.006838000957720374, 0.012843416735794254, 0.011257735386877264}};

    // heavy tails, fraction 1: passes 8 converged 0, core 95504 far 4496 halo 2748
    const ExpectedSelection heavy100Expected =
        {8, false, 100000, 95504, 4496, 2748, 1.0000000000000003e-10, 4.496000000000002e-12,
         15719.20857347983, {8.763062603691254e-06, -7.804185895869224e-06, -5.741339299266207e-06},
         {0.0023586196350638086, 0.002999543015084024, 0.0025019119017746644},
         {-0.014136092688360932, -0.012501022715932557, -0.010904547591835874},
         {0.014157442373975399, 0.013903888817967218, 0.011257735386877264}};

    // scraped, nSigma 12: passes 3 converged 1, core 201979 far 21 halo 95
    const ExpectedSelection scraped12Expected =
        {3, true, 202000, 201979, 21, 95, 2.0199999999999997e-10, 2.0999999999999996e-14,
         14.560080083496405,
         {-0.00013782758998007722, -6.16470676303518e-06, 1.4220842864287944e-06},
         {0.0008847595290983996, 0.0030115993205686808, 0.0025127535846909146},
         {-0.010680040906320612, -0.035609541878147126, -0.02993666639744278},
         {0.010226145866210179, 0.03571405145445878, 0.027596287644470228}};

    // the Gaussian, flat at z = 0: passes 2 converged 1, core 200001 far 39 halo 36
    const ExpectedSelection flatExpected =
        {2, true, 200040, 200001, 39, 36, 2.0003999999999996e-10, 3.9000000000000017e-14,
         34.89978169569404, {3.543148269927127e-07, 1.5885330220825716e-06, 0.0},
         {0.000995974925626115, 0.0030005915692836757, 0.0},
         {-0.004730037975479777, -0.013018094678226786, 0.0},
         {0.004442721171545147, 0.013430044993129886, 0.0}};

    // the Gaussian, core flat at z = 0: passes 2 converged 1, core 200000 far 40 halo 40
    const ExpectedSelection flatCoreExpected =
        {2, true, 200040, 200000, 40, 40, 2.0003999999999996e-10, 4.0000000000000006e-14,
         34.89970691537061, {3.560613175506539e-07, 1.5326512139941648e-06, 0.0},
         {0.0009959771097653059, 0.0030004950150713705, 0.0},
         {-0.004730037975479777, -0.013018094678226786, 0.0},
         {0.004442721171545147, 0.013430044993129886, 0.0}};

    // the Gaussian without halo, flat at z = 0: passes 1 converged 1, core 200000 far 0 halo 0
    const ExpectedSelection flatPureExpected =
        {1, true, 200000, 200000, 0, 0, 1.9999999999999996e-10, 0.0, 0.0,
         {3.520608387473155e-07, 1.4726061329464281e-06, 0.0},
         {0.0009959762260057924, 0.0030006076512576426, 0.0},
         {-0.004730037975479777, -0.013018094678226786, 0.0},
         {0.004442721171545147, 0.013430044993129886, 0.0}};

    const double pythonTolerance = 1e-10;

    Vector_t toVector(const double v[3]) {
        return Vector_t({v[0], v[1], v[2]});
    }

    /// |actual - expected| <= tolerance * |expected|
    void expectVectorNear(const Vector_t& actual, const Vector_t& expected, double tolerance,
                          const std::string& what) {
        EXPECT_LE(euclidean_norm(actual - expected), tolerance * euclidean_norm(expected))
            << what << ": " << actual << " against " << expected;
    }

    void expectRelativeNear(double actual, double expected, double tolerance,
                            const std::string& what) {
        EXPECT_LE(std::abs(actual - expected), tolerance * std::abs(expected))
            << what << ": " << actual << " against " << expected;
    }

    /// The C++ selection against the Python reference
    void compareSelectionWithPython(const CoreSelection& s, const ExpectedSelection& e) {
        const double tol = pythonTolerance;
        EXPECT_EQ(s.passes, e.passes);
        EXPECT_EQ(s.converged, e.converged);
        EXPECT_EQ(s.numLive, e.numLive);
        EXPECT_EQ(s.numCore, e.numCore);
        EXPECT_EQ(s.numFar, e.numFar);
        EXPECT_EQ(s.numHalo, e.numHalo);
        expectRelativeNear(s.chargeLive, e.chargeLive, tol, "chargeLive");
        expectRelativeNear(s.chargeFar, e.chargeFar, tol, "chargeFar");
        expectRelativeNear(s.farMaxDistance, e.farMaxDistance, tol, "farMaxDistance");
        for (unsigned int d = 0; d < 3; ++d) {
            // the mean is close to 0: compare it, and sigma, on the scale of sigma (at least
            // flatWidth, for a flat axis)
            const double scale = std::max(e.sigma[d], flatWidth);
            EXPECT_LE(std::abs(s.mean[d] - e.mean[d]), tol * scale) << "mean " << d;
            EXPECT_LE(std::abs(s.sigma[d] - e.sigma[d]), tol * scale) << "sigma " << d;
            // the bounds are particle coordinates: exact
            EXPECT_EQ(s.boundsMin[d], e.boundsMin[d]) << "boundsMin " << d;
            EXPECT_EQ(s.boundsMax[d], e.boundsMax[d]) << "boundsMax " << d;
        }
    }

    /// The C++ selection and far-field model against the Python reference
    void compareWithPython(const Bunch& b, const Expected& e) {
        const FarModel m = computeFarModel(b, syntheticGamma);
        const CoreSelection& s = m.selection;
        const double tol = pythonTolerance;
        compareSelectionWithPython(s, e.selection);

        // classify() puts exactly the selected core into CORE; the far table holds the same
        // particles in the same classes, and the same rows are flagged
        const std::size_t numCoreClass = std::count(m.classes.begin(), m.classes.end(),
                                                    ParticleClass::CORE);
        EXPECT_EQ(numCoreClass, s.numCore);
        EXPECT_EQ(m.table.size(), e.tableSize);
        EXPECT_EQ(m.numMesh, e.numMesh);
        std::vector<std::size_t> codes;
        for (std::size_t i: m.table) {
            codes.push_back(4 * i + static_cast<std::size_t>(m.classes[i]));
        }
        EXPECT_EQ(hashIndices(codes), e.tableHash);
        EXPECT_EQ(m.flagged.size(), e.numFlagged);
        EXPECT_EQ(hashIndices(m.flagged), e.flaggedHash);

        expectRelativeNear(m.moments.charge, e.coreCharge, tol, "core charge");
        // the centroid is close to the centre: compare it on the scale of the rest-frame sigma
        const Vector_t restSigma({s.sigma[0], syntheticGamma * s.sigma[1], s.sigma[2]});
        EXPECT_LE(euclidean_norm(m.moments.centroid - toVector(e.centroid)),
                  tol * euclidean_norm(restSigma))
            << "centroid " << m.moments.centroid;
        const double qNorm = std::max({std::abs(e.quadrupole[0]), std::abs(e.quadrupole[3]),
                                       std::abs(e.quadrupole[5])});
        unsigned int l = 0;
        for (unsigned int i = 0; i < 3; ++i) {
            for (unsigned int j = i; j < 3; ++j, ++l) {
                EXPECT_LE(std::abs(m.moments.quadrupole[i][j] - e.quadrupole[l]), tol * qNorm)
                    << "quadrupole " << i << j;
                EXPECT_EQ(m.moments.quadrupole[i][j], m.moments.quadrupole[j][i]);
            }
        }
        const double trace = (m.moments.quadrupole[0][0] + m.moments.quadrupole[1][1] +
                              m.moments.quadrupole[2][2]);
        EXPECT_LE(std::abs(trace), 1e-14 * qNorm);
        expectRelativeNear(m.eps, e.eps, tol, "eps");

        std::vector<Vector_t> targets;
        for (unsigned int k = 0; k < 4; ++k) {
            const std::size_t row = e.rows[k];
            ASSERT_LT(row, m.table.size());
            expectVectorNear(m.multipole[row], toVector(e.multipole[k]), tol, "multipole");
            expectRelativeNear(m.indicator[row], e.indicator[k], tol, "indicator");
            expectVectorNear(m.farFar[row], toVector(e.farFar[k]), tol, "far-far");
            targets.push_back(m.position[row]);
        }
        std::vector<Vector_t> coreField(4, Vector_t(0.0));
        addDirectField(targets, 0, 4, m.coreSource, m.coreCharge, 0.0, coreField);
        for (unsigned int k = 0; k < 4; ++k) {
            expectVectorNear(coreField[k], toVector(e.coreField[k]), tol, "core field");
        }

        double sums[4] = {0.0, 0.0, 0.0, 0.0};
        for (std::size_t i = 0; i < m.table.size(); ++i) {
            sums[0] += euclidean_norm(m.multipole[i]);
            sums[1] += euclidean_norm(m.farFar[i]);
            sums[2] += m.indicator[i];
            sums[3] = std::max(sums[3], m.indicator[i]);
        }
        for (unsigned int k = 0; k < 4; ++k) {
            expectRelativeNear(sums[k], e.sums[k], tol, "table sum");
        }
    }

    /// MAD-equivalent sigma of one coordinate
    double madSigma(std::vector<double> x) {
        auto median = [](std::vector<double>& v) {
            const std::size_t n = v.size();
            std::nth_element(v.begin(), v.begin() + n / 2, v.end());
            const double upper = v[n / 2];
            if (n % 2 == 1) {
                return upper;
            }
            return 0.5 * (upper + *std::max_element(v.begin(), v.begin() + n / 2));
        };
        const double m = median(x);
        for (double& v: x) {
            v = std::abs(v - m);
        }
        return 1.482602218505602 * median(x);
    }

    /// Simulated ranks for the tests: every thread is a rank, and each reduction combines the
    /// ranks' values in rank order, so every rank gets the same result.
    class ThreadReducer: public Reducer {
    public:
        struct Shared {
            explicit Shared(int n): numRanks(n), barrier(n), slots(n) { }
            int numRanks;
            std::barrier<> barrier;
            std::vector<std::vector<double>> slots;
            std::atomic<bool> mismatch{false};
            std::atomic<unsigned int> numCalls{0};
        };

        ThreadReducer(Shared& shared, int rank): shared_m(shared), rank_m(rank) { }

        void sum(double* values, unsigned int count) const override {
            combine(values, count, [](double a, double b) { return a + b; });
        }

        void min(double* values, unsigned int count) const override {
            combine(values, count, [](double a, double b) { return std::min(a, b); });
        }

    private:
        void combine(double* values, unsigned int count,
                     const std::function<double(double, double)>& op) const {
            shared_m.slots[rank_m].assign(values, values + count);
            if (rank_m == 0) {
                ++ shared_m.numCalls;
            }
            shared_m.barrier.arrive_and_wait();
            for (int p = 0; p < shared_m.numRanks; ++p) {
                if (shared_m.slots[p].size() != count) {
                    shared_m.mismatch = true;
                    count = std::min<unsigned int>(count, shared_m.slots[p].size());
                }
            }
            for (unsigned int k = 0; k < count; ++k) {
                double v = shared_m.slots[0][k];
                for (int p = 1; p < shared_m.numRanks; ++p) {
                    v = op(v, shared_m.slots[p][k]);
                }
                values[k] = v;
            }
            shared_m.barrier.arrive_and_wait();
        }

        Shared& shared_m;
        int rank_m;
    };

    /// Runs body(rank, reducer) on numRanks simulated ranks; returns the number of
    /// reductions, which every rank must have made with the same counts
    unsigned int runOnRanks(int numRanks,
                            const std::function<void(int, const Reducer&)>& body) {
        ThreadReducer::Shared shared(numRanks);
        std::vector<std::thread> ranks;
        for (int p = 0; p < numRanks; ++p) {
            ranks.emplace_back([&shared, &body, p]() {
                const ThreadReducer reducer(shared, p);
                body(p, reducer);
            });
        }
        for (std::thread& t: ranks) {
            t.join();
        }
        EXPECT_FALSE(shared.mismatch);
        return shared.numCalls.load();
    }

    /// Splits a bunch over numRanks: rank 0 gets every particle outside 5 sigma of the
    /// Gaussian (the halo, the tails and the lost ones) and a small share of the rest, the
    /// others unequal shares.
    std::vector<Bunch> distribute(const Bunch& b, int numRanks) {
        std::vector<Bunch> parts(numRanks);
        std::vector<std::size_t> inner;
        for (std::size_t i = 0; i < b.R.size(); ++i) {
            bool outer = false;
            for (unsigned int d = 0; d < 3; ++d) {
                if (std::abs(b.R[i][d]) > 5.0 * syntheticSigma[d]) {
                    outer = true;
                }
            }
            if (outer) {
                parts[0].push(b.R[i], b.Q[i], b.bin[i]);
            } else {
                inner.push_back(i);
            }
        }
        // weights 1, 2, ..., numRanks
        const std::size_t total = numRanks * (numRanks + 1) / 2;
        std::size_t start = 0;
        for (int p = 0; p < numRanks; ++p) {
            const std::size_t end = (p + 1 == numRanks ? inner.size() :
                                     start + inner.size() * (p + 1) / total);
            for (std::size_t k = start; k < end; ++k) {
                const std::size_t i = inner[k];
                parts[p].push(b.R[i], b.Q[i], b.bin[i]);
            }
            start = end;
        }
        return parts;
    }

    /// Particle i on rank i % numRanks
    std::vector<Bunch> distributeRoundRobin(const Bunch& b, int numRanks) {
        std::vector<Bunch> parts(numRanks);
        for (std::size_t i = 0; i < b.R.size(); ++i) {
            parts[i % numRanks].push(b.R[i], b.Q[i], b.bin[i]);
        }
        return parts;
    }

    void expectSameSelection(const CoreSelection& a, const CoreSelection& b, double tolerance) {
        EXPECT_EQ(a.passes, b.passes);
        EXPECT_EQ(a.converged, b.converged);
        EXPECT_EQ(a.numLive, b.numLive);
        EXPECT_EQ(a.numCore, b.numCore);
        EXPECT_EQ(a.numFar, b.numFar);
        EXPECT_EQ(a.numHalo, b.numHalo);
        expectRelativeNear(a.chargeLive, b.chargeLive, tolerance, "chargeLive");
        expectRelativeNear(a.chargeFar, b.chargeFar, tolerance, "chargeFar");
        expectRelativeNear(a.farMaxDistance, b.farMaxDistance, tolerance, "farMaxDistance");
        for (unsigned int d = 0; d < 3; ++d) {
            EXPECT_LE(std::abs(a.mean[d] - b.mean[d]), tolerance * b.sigma[d]);
            expectRelativeNear(a.sigma[d], b.sigma[d], tolerance, "sigma");
            EXPECT_EQ(a.boundsMin[d], b.boundsMin[d]);
            EXPECT_EQ(a.boundsMax[d], b.boundsMax[d]);
        }
    }

    void expectSameMoments(const CoreMoments& a, const CoreMoments& b, double tolerance) {
        expectRelativeNear(a.charge, b.charge, tolerance, "charge");
        expectVectorNear(a.centroid, b.centroid, tolerance, "centroid");
        for (unsigned int i = 0; i < 3; ++i) {
            expectVectorNear(a.quadrupole[i], b.quadrupole[i], tolerance, "quadrupole");
        }
    }

    /// Exact field of the point charges at t in long double, without the coupling constant
    Vektor<long double, 3> exactField(const Vector_t& t, const std::vector<Vector_t>& s,
                                      const std::vector<double>& q) {
        Vektor<long double, 3> E(0.0L);
        for (std::size_t j = 0; j < s.size(); ++j) {
            Vektor<long double, 3> d;
            for (unsigned int k = 0; k < 3; ++k) {
                d[k] = static_cast<long double>(t[k]) - static_cast<long double>(s[j][k]);
            }
            const long double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            const long double w = q[j] / (r2 * std::sqrt(r2));
            for (unsigned int k = 0; k < 3; ++k) {
                E[k] += w * d[k];
            }
        }
        return E;
    }

    /// Target directions: the axes, the face and space diagonals and a few odd ones
    std::vector<Vector_t> directions() {
        std::vector<Vector_t> dirs;
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    if (i != 0 || j != 0 || k != 0) {
                        const Vector_t v({double(i), double(j), double(k)});
                        dirs.push_back(v / euclidean_norm(v));
                    }
                }
            }
        }
        for (const Vector_t& v: {Vector_t({0.3, -0.7, 0.2}), Vector_t({-0.9, 0.1, 0.45}),
                                 Vector_t({0.05, 0.2, -1.0})}) {
            dirs.push_back(v / euclidean_norm(v));
        }
        return dirs;
    }
}


// Selection on the Gaussian with 40 halo particles
TEST(CoreFitSCTest, SelectionGaussian) {
    const Bunch b = makeSynthetic(false);
    const CoreSelection s = select(b);

    EXPECT_EQ(s.numLive, syntheticCore + 40);     // the 10 lost particles count nowhere
    EXPECT_EQ(s.numFar, 40u);
    EXPECT_EQ(s.numHalo, 40u);
    EXPECT_EQ(s.numCore, syntheticCore);
    EXPECT_TRUE(s.converged);
    EXPECT_LE(s.passes, 6u);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_LT(std::abs(s.sigma[d] / syntheticSigma[d] - 1.0), 5e-3) << "axis " << d;
    }
    expectRelativeNear(s.getFarChargeFraction(), 40.0 / (syntheticCore + 40), 1e-10,
                       "far charge fraction");
    EXPECT_GT(s.farMaxDistance, 15.0);
    EXPECT_LT(s.farMaxDistance, 35.0);

    compareWithPython(b, gaussianExpected);
}


// A scraped bunch (Gaussian truncated at +1.5 sigma in x, plus 1 % tails)
TEST(CoreFitSCTest, SelectionScraped) {
    const Bunch b = makeSynthetic(true);
    const CoreSelection s = select(b);

    // the core part: the truncated Gaussian without the tails
    std::vector<double> x[3];
    double rms[3], rmsLive[3];
    for (unsigned int d = 0; d < 3; ++d) {
        double s1 = 0.0, s2 = 0.0, t1 = 0.0, t2 = 0.0, nLive = 0.0;
        for (std::size_t i = 0; i < b.R.size(); ++i) {
            if (i < syntheticCore) {
                x[d].push_back(b.R[i][d]);
                s1 += b.R[i][d];
                s2 += b.R[i][d] * b.R[i][d];
            }
            if (b.bin[i] >= 0) {
                t1 += b.R[i][d];
                t2 += b.R[i][d] * b.R[i][d];
                nLive += 1.0;
            }
        }
        rms[d] = std::sqrt(s2 / syntheticCore - std::pow(s1 / syntheticCore, 2));
        rmsLive[d] = std::sqrt(t2 / nLive - std::pow(t1 / nLive, 2));
    }
    for (unsigned int d = 0; d < 3; ++d) {
        const double mad = madSigma(x[d]);
        // within 5 % of the MAD-equivalent sigma of the core part (x: 4.5 %; the rms of a
        // Gaussian cut at +1.5 sigma is itself 5.1 % below its MAD sigma)
        EXPECT_LT(std::abs(s.sigma[d] / mad - 1.0), 0.05) << "axis " << d;
        // the tails inflate the plain rms by 7-9 %, the clipped sigma stays within 1 % of the
        // rms of the core part
        EXPECT_GT(rmsLive[d] / rms[d], 1.05) << "axis " << d;
        EXPECT_LT(std::abs(s.sigma[d] / rms[d] - 1.0), 0.01) << "axis " << d;
    }

    compareWithPython(b, scrapedExpected);
}


// Heavy tails, x a Cauchy variable in 25 % and in all of the particles. The clipped passes
// shrink sigma_x slowly: the first bunch converges on the last pass, the second stops there
// unconverged, as the selection does at turns 119 and 120 of the IsoDAR run.
TEST(CoreFitSCTest, SelectionHeavyTails) {
    const CoreSelection s25 = select(makeHeavyTailed(0.25));
    EXPECT_EQ(s25.passes, maxClipPasses);
    EXPECT_TRUE(s25.converged);
    compareSelectionWithPython(s25, heavy25Expected);

    const CoreSelection s100 = select(makeHeavyTailed(1.0));
    EXPECT_EQ(s100.passes, maxClipPasses);
    EXPECT_FALSE(s100.converged);
    compareSelectionWithPython(s100, heavy100Expected);
}


// nSigma only sets the final pass. With nSigma = 12 the halo still counts every live
// particle beyond 10 sigma, core particles included; with nSigma = 1e6 no particle is far
// and the bounds are the extremes of the live particles.
TEST(CoreFitSCTest, SelectionCoreWidth) {
    const Bunch b = makeSynthetic(true);
    const CoreSelection s6 = select(b);

    const CoreSelection s12 = select(b, SerialReducer(), 12.0);
    compareSelectionWithPython(s12, scraped12Expected);
    EXPECT_EQ(s12.passes, s6.passes);
    EXPECT_EQ(s12.numHalo, s6.numHalo);
    EXPECT_GT(s12.numHalo, s12.numFar);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_EQ(s12.sigma[d], s6.sigma[d]);
    }

    const CoreSelection all = select(b, SerialReducer(), 1e6);
    EXPECT_EQ(all.numFar, 0u);
    EXPECT_EQ(all.numCore, all.numLive);
    EXPECT_EQ(all.numHalo, s6.numHalo);
    EXPECT_EQ(all.chargeFar, 0.0);
    EXPECT_EQ(all.farMaxDistance, 0.0);
    Vector_t lower(std::numeric_limits<double>::max()), upper(-std::numeric_limits<double>::max());
    for (std::size_t i = 0; i < b.R.size(); ++i) {
        if (b.bin[i] >= 0) {
            for (unsigned int d = 0; d < 3; ++d) {
                lower[d] = std::min(lower[d], b.R[i][d]);
                upper[d] = std::max(upper[d], b.R[i][d]);
            }
        }
    }
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_EQ(all.boundsMin[d], lower[d]);
        EXPECT_EQ(all.boundsMax[d], upper[d]);
    }
}


TEST(CoreFitSCTest, SelectionDegenerate) {
    // no particle
    CoreSelection s = selectCore(nullptr, nullptr, nullptr, 0, 6.0, 4.0, SerialReducer());
    EXPECT_EQ(s.numLive, 0u);
    EXPECT_EQ(s.numCore, 0u);
    EXPECT_EQ(s.numFar, 0u);
    EXPECT_EQ(s.passes, 0u);
    EXPECT_EQ(s.getFarChargeFraction(), 0.0);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_GT(s.boundsMin[d], s.boundsMax[d]);
        EXPECT_EQ(s.mean[d], 0.0);
        EXPECT_EQ(s.sigma[d], 0.0);
    }
    EXPECT_FALSE(isInside(Vector_t(0.0), s.boundsMin, s.boundsMax));

    // only lost particles
    Bunch b;
    b.R = {Vector_t({1e-3, 2e-3, 3e-3}), Vector_t({-1e-3, 0.0, 1e-3})};
    b.Q = {1e-15, 1e-15};
    b.bin = {-1, -1};
    s = select(b);
    EXPECT_EQ(s.numLive, 0u);
    EXPECT_EQ(s.numCore, 0u);
    EXPECT_EQ(s.chargeLive, 0.0);

    // one live particle: it is the core
    b.bin = {0, -1};
    s = select(b);
    EXPECT_EQ(s.numLive, 1u);
    EXPECT_EQ(s.numCore, 1u);
    EXPECT_EQ(s.numFar, 0u);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_EQ(s.boundsMin[d], b.R[0][d]);
        EXPECT_EQ(s.boundsMax[d], b.R[0][d]);
    }

    // two particles exactly one sigma from the mean in every axis (powers of 2, so the sums
    // are exact): the clip and core tests include their boundary
    const double a = std::ldexp(1.0, -10);
    b.R = {Vector_t({a, 2.0 * a, 4.0 * a}), Vector_t({-a, -2.0 * a, -4.0 * a})};
    b.bin = {0, 0};
    s = selectCore(b.R.data(), b.Q.data(), b.bin.data(), 2, 1.0, 1.0, SerialReducer());
    EXPECT_EQ(s.passes, 1u);
    EXPECT_TRUE(s.converged);
    EXPECT_EQ(s.numCore, 2u);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_EQ(s.sigma[d], b.R[0][d]);
        EXPECT_EQ(s.boundsMax[d], b.R[0][d]);
    }
    // a clip at half a sigma selects nothing: that pass is dropped, the plain rms stays
    s = selectCore(b.R.data(), b.Q.data(), b.bin.data(), 2, 1.0, 0.5, SerialReducer());
    EXPECT_EQ(s.passes, 0u);
    EXPECT_FALSE(s.converged);
    EXPECT_EQ(s.numCore, 2u);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_EQ(s.sigma[d], b.R[0][d]);
    }
}


// A bunch flat in z (every z the same, as in a bunch tracked in the median plane), at z = 0
// and at offsets where the mean of z is not exactly z0. With the smallest half-width
// flatWidth, the round-off of the mean does not decide whether the plane is in the core, so
// every offset and rank count gives the selection of the other two axes. With only the core
// in the plane, the particles off it are far; their distance counts in x and y only. The
// Gaussian alone, without the halo, changes its count on the first clipped pass, so there
// the stop test decides the pass count: a flat axis has settled by its new rms alone, since
// the rms of pass 0, summed about the origin, is round-off that depends on the offset and
// on the order of the sum.
TEST(CoreFitSCTest, SelectionFlatAxis) {
    const Bunch gaussian = makeSynthetic(false);
    // the Gaussian particles alone, without the halo and the lost particles
    Bunch pure = gaussian;
    pure.R.resize(syntheticCore);
    pure.Q.resize(syntheticCore);
    pure.bin.resize(syntheticCore);
    const struct {
        const Bunch* bunch;
        bool coreOnly;
        const ExpectedSelection* expected;
        std::string label;
    } cases[] = {{&gaussian, false, &flatExpected, "all flat"},
                 {&gaussian, true, &flatCoreExpected, "core flat"},
                 {&pure, false, &flatPureExpected, "no halo, all flat"}};
    for (const auto& c: cases) {
        const CoreSelection s0 = select(makeFlat(*c.bunch, 0.0, c.coreOnly));
        compareSelectionWithPython(s0, *c.expected);
        if (c.coreOnly) {
            EXPECT_EQ(s0.numCore, syntheticCore);
            EXPECT_EQ(s0.numFar, 40u);
            EXPECT_EQ(s0.numHalo, 40u);
        }
        for (double z0: {0.0, 1e-3, 0.1, 1.9}) {
            const Bunch b = makeFlat(*c.bunch, z0, c.coreOnly);
            for (int numRanks: {1, 3, 8}) {
                const std::vector<Bunch> parts = distributeRoundRobin(b, numRanks);
                std::vector<CoreSelection> result(numRanks);
                runOnRanks(numRanks, [&](int p, const Reducer& reducer) {
                    result[p] = select(parts[p], reducer);
                });
                for (const CoreSelection& s: result) {
                    SCOPED_TRACE(c.label + ", z0 " + std::to_string(z0) + ", " +
                                 std::to_string(numRanks) + " ranks");
                    EXPECT_EQ(s.passes, s0.passes);
                    EXPECT_EQ(s.converged, s0.converged);
                    EXPECT_EQ(s.numCore, s0.numCore);
                    EXPECT_EQ(s.numFar, s0.numFar);
                    EXPECT_EQ(s.numHalo, s0.numHalo);
                    expectRelativeNear(s.farMaxDistance, s0.farMaxDistance, 1e-10,
                                       "farMaxDistance");
                    for (unsigned int d = 0; d < 2; ++d) {
                        expectRelativeNear(s.sigma[d], s0.sigma[d], 1e-10, "sigma");
                        EXPECT_EQ(s.boundsMin[d], s0.boundsMin[d]);
                        EXPECT_EQ(s.boundsMax[d], s0.boundsMax[d]);
                    }
                    EXPECT_LE(s.sigma[2], flatWidth);
                    EXPECT_EQ(s.boundsMin[2], z0);
                    EXPECT_EQ(s.boundsMax[2], z0);
                }
            }
        }
    }

    // a core in a slab of +-1e-12 m, with an rms above 0 but below flatWidth: the same
    const Bunch flatCore = makeFlat(gaussian, 0.0, true);
    Bunch slab = flatCore;
    for (std::size_t i = 0; i < syntheticCore; ++i) {
        slab.R[i][2] = (i % 2 == 0 ? 1e-12 : -1e-12);
    }
    const CoreSelection sFlat = select(flatCore);
    const CoreSelection sSlab = select(slab);
    EXPECT_GT(sSlab.sigma[2], 0.0);
    EXPECT_LE(sSlab.sigma[2], flatWidth);
    EXPECT_EQ(sSlab.numCore, sFlat.numCore);
    EXPECT_EQ(sSlab.numFar, sFlat.numFar);
    EXPECT_EQ(sSlab.numHalo, sFlat.numHalo);
    expectRelativeNear(sSlab.farMaxDistance, sFlat.farMaxDistance, 1e-10, "farMaxDistance");
    EXPECT_EQ(sSlab.boundsMin[2], -1e-12);
    EXPECT_EQ(sSlab.boundsMax[2], 1e-12);
}


// The reduction order: P = 1, 3 and 8 simulated ranks, with every particle outside 5 sigma
// on rank 0, give the same selection up to round-off, the same counts and exact bounds.
TEST(CoreFitSCTest, SelectionIndependentOfRankCount) {
    const std::vector<Bunch> bunches = {makeSynthetic(false), makeSynthetic(true),
                                        makeHeavyTailed(1.0)};
    for (const Bunch& b: bunches) {
        const CoreSelection serial = select(b);
        const CoreMoments serialMoments = getCoreMoments(b, serial, SerialReducer());

        for (int numRanks: {1, 3, 8}) {
            const std::vector<Bunch> parts = distribute(b, numRanks);
            std::vector<CoreSelection> result(numRanks);
            std::vector<CoreMoments> moments(numRanks);
            const unsigned int numCalls = runOnRanks(numRanks, [&](int p, const Reducer& reducer) {
                result[p] = select(parts[p], reducer);
                moments[p] = getCoreMoments(parts[p], result[p], reducer);
            });
            // pass 0, the clipped passes, the far sums and the extrema, then the moments
            EXPECT_EQ(numCalls, serial.passes + 4);
            // round-off: sums of 2e5 terms in a different order differ by up to N eps = 2e-11
            for (int p = 0; p < numRanks; ++p) {
                expectSameSelection(result[p], serial, 1e-10);
                expectSameMoments(moments[p], serialMoments, 1e-10);
                // every rank gets bitwise the same result
                EXPECT_EQ(result[p].sigma[0], result[0].sigma[0]);
                EXPECT_EQ(moments[p].quadrupole[2][2], moments[0].quadrupole[2][2]);
            }
        }
    }
}


// The same with IpplReducer on all ranks of the test. On one rank the reductions do nothing:
// run with mpirun -np 4 (or more) for a real test.
TEST(CoreFitSCTest, SelectionIpplReducer) {
    for (bool scraped: {false, true}) {
        const Bunch b = makeSynthetic(scraped);
        const CoreSelection serial = select(b);
        const CoreMoments serialMoments = getCoreMoments(b, serial, SerialReducer());

        const std::size_t numNodes = Ippl::getNodes();
        const std::size_t node = Ippl::myNode();
        const std::size_t n = b.R.size();
        Bunch part;
        for (std::size_t i = node * n / numNodes; i < (node + 1) * n / numNodes; ++i) {
            part.push(b.R[i], b.Q[i], b.bin[i]);
        }
        const CoreSelection s = select(part, IpplReducer());
        expectSameSelection(s, serial, 1e-10);
        expectSameMoments(getCoreMoments(part, s, IpplReducer()), serialMoments, 1e-10);
    }
}


// The box of the core-fitted mesh: the core bounds, widened to the particles flagged lost
// in the core window, and never by those outside it. With nSigma 1e6 it is the bounds of all
// particles, lost ones included, as MESHFIT="ALL" fits the mesh. A live particle lies in the
// box exactly when it is a core particle, and the box does not depend on the rank count.
TEST(CoreFitSCTest, MeshBox) {
    // the 10 lost particles of the synthetic bunch lie at 20-50 sigma, outside the window
    Bunch b = makeSynthetic(false);
    const CoreSelection s = select(b);
    Vector_t boxMin, boxMax;
    getMeshBox(b.R.data(), b.bin.data(), b.R.size(), s, 6.0, SerialReducer(), boxMin, boxMax);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_EQ(boxMin[d], s.boundsMin[d]);
        EXPECT_EQ(boxMax[d], s.boundsMax[d]);
    }

    // lost particles in the window beyond the core bounds (at most 5 sigma here) widen the
    // box; those outside the window in one axis do not, nor do they change the selection
    const Vector_t& mean = s.mean;
    const Vector_t& sigma = s.sigma;
    for (unsigned int d = 0; d < 3; ++d) {
        ASSERT_LT(s.boundsMax[d], mean[d] + 5.0 * sigma[d]);
        ASSERT_GT(s.boundsMin[d], mean[d] - 5.0 * sigma[d]);
    }
    const Vector_t widenX({mean[0] + 5.9 * sigma[0], mean[1], mean[2]});
    const Vector_t widenYZ({mean[0], mean[1] - 5.9 * sigma[1], mean[2] + 5.9 * sigma[2]});
    b.push(widenX, syntheticCharge, -1);
    b.push(widenYZ, syntheticCharge, -1);
    b.push(Vector_t({mean[0] - 6.1 * sigma[0], mean[1], mean[2]}), syntheticCharge, -1);
    b.push(Vector_t({mean[0], mean[1] + 5.9 * sigma[1], mean[2] - 6.1 * sigma[2]}),
           syntheticCharge, -1);
    const CoreSelection sLost = select(b);
    expectSameSelection(sLost, s, 0.0);
    getMeshBox(b.R.data(), b.bin.data(), b.R.size(), sLost, 6.0, SerialReducer(),
               boxMin, boxMax);
    EXPECT_EQ(boxMin[0], s.boundsMin[0]);
    EXPECT_EQ(boxMax[0], widenX[0]);
    EXPECT_EQ(boxMin[1], widenYZ[1]);
    EXPECT_EQ(boxMax[1], s.boundsMax[1]);
    EXPECT_EQ(boxMin[2], s.boundsMin[2]);
    EXPECT_EQ(boxMax[2], widenYZ[2]);

    // classes on the box: CORE exactly for the live particles in the window, INSIDE for the
    // two lost ones in it
    std::size_t numCore = 0, numInside = 0;
    for (std::size_t i = 0; i < b.R.size(); ++i) {
        bool inWindow = true;
        for (unsigned int d = 0; d < 3; ++d) {
            if (std::abs(b.R[i][d] - mean[d]) > 6.0 * sigma[d]) {
                inWindow = false;
            }
        }
        const ParticleClass c = classify(b.bin[i], b.R[i], boxMin, boxMax, boxMin, boxMax);
        if (c == ParticleClass::CORE) {
            ++ numCore;
        } else if (c == ParticleClass::INSIDE) {
            ++ numInside;
        }
        const ParticleClass inBox = (b.bin[i] >= 0 ? ParticleClass::CORE :
                                     ParticleClass::INSIDE);
        EXPECT_EQ(c == inBox, inWindow) << "particle " << i;
    }
    EXPECT_EQ(numCore, s.numCore);
    EXPECT_EQ(numInside, 2u);

    // nSigma 1e6: the bounds of all particles, live and lost
    const CoreSelection sAll = select(b, SerialReducer(), 1e6);
    getMeshBox(b.R.data(), b.bin.data(), b.R.size(), sAll, 1e6, SerialReducer(),
               boxMin, boxMax);
    Vector_t lower(std::numeric_limits<double>::max()), upper(-std::numeric_limits<double>::max());
    for (const Vector_t& r: b.R) {
        for (unsigned int d = 0; d < 3; ++d) {
            lower[d] = std::min(lower[d], r[d]);
            upper[d] = std::max(upper[d], r[d]);
        }
    }
    bool lostOutsideLive = false;
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_EQ(boxMin[d], lower[d]);
        EXPECT_EQ(boxMax[d], upper[d]);
        if (lower[d] < sAll.boundsMin[d] || upper[d] > sAll.boundsMax[d]) {
            lostOutsideLive = true;
        }
    }
    // (the lost particles at 20-50 sigma lie beyond the live ones in some axis)
    EXPECT_TRUE(lostOutsideLive);

    // 1, 3 and 8 simulated ranks, the lost particles spread over them: the same box
    getMeshBox(b.R.data(), b.bin.data(), b.R.size(), sLost, 6.0, SerialReducer(),
               boxMin, boxMax);
    for (int numRanks: {1, 3, 8}) {
        const std::vector<Bunch> parts = distributeRoundRobin(b, numRanks);
        std::vector<Vector_t> rankMin(numRanks), rankMax(numRanks);
        const unsigned int numCalls = runOnRanks(numRanks, [&](int p, const Reducer& reducer) {
            const Bunch& part = parts[p];
            const CoreSelection sp = select(part, reducer);
            getMeshBox(part.R.data(), part.bin.data(), part.R.size(), sp, 6.0, reducer,
                       rankMin[p], rankMax[p]);
        });
        // the selection, then one minimum
        EXPECT_EQ(numCalls, sLost.passes + 4);
        for (int p = 0; p < numRanks; ++p) {
            for (unsigned int d = 0; d < 3; ++d) {
                EXPECT_EQ(rankMin[p][d], boxMin[d]);
                EXPECT_EQ(rankMax[p][d], boxMax[d]);
            }
        }
    }

    // and on all ranks of the test, through IpplReducer (mpirun -np 4 for a real test)
    const std::size_t numNodes = Ippl::getNodes();
    const std::size_t node = Ippl::myNode();
    Bunch part;
    for (std::size_t i = node; i < b.R.size(); i += numNodes) {
        part.push(b.R[i], b.Q[i], b.bin[i]);
    }
    const CoreSelection sPart = select(part, IpplReducer());
    Vector_t partMin, partMax;
    getMeshBox(part.R.data(), part.bin.data(), part.R.size(), sPart, 6.0, IpplReducer(),
               partMin, partMax);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_EQ(partMin[d], boxMin[d]);
        EXPECT_EQ(partMax[d], boxMax[d]);
    }

    // the window includes its boundary: two live particles one sigma from the mean (powers of
    // 2, so the mean and rms are exact), nSigma 2
    const double a = std::ldexp(1.0, -10);
    Bunch pair;
    pair.push(Vector_t({a, 2.0 * a, 4.0 * a}), 1e-15, 0);
    pair.push(Vector_t({-a, -2.0 * a, -4.0 * a}), 1e-15, 0);
    pair.push(Vector_t({2.0 * a, 0.0, 0.0}), 1e-15, -1);
    pair.push(Vector_t({0.0, -std::nextafter(4.0 * a, 1.0), 0.0}), 1e-15, -1);
    const CoreSelection sPair = selectCore(pair.R.data(), pair.Q.data(), pair.bin.data(), 4,
                                           2.0, 1.0, SerialReducer());
    ASSERT_EQ(sPair.sigma[0], a);
    getMeshBox(pair.R.data(), pair.bin.data(), 4, sPair, 2.0, SerialReducer(),
               boxMin, boxMax);
    EXPECT_EQ(boxMax[0], 2.0 * a);
    EXPECT_EQ(boxMin[1], -2.0 * a);

    // no particle: an empty box
    const CoreSelection sNone = selectCore(nullptr, nullptr, nullptr, 0, 6.0, 4.0,
                                           SerialReducer());
    getMeshBox(nullptr, nullptr, 0, sNone, 6.0, SerialReducer(), boxMin, boxMax);
    for (unsigned int d = 0; d < 3; ++d) {
        EXPECT_GT(boxMin[d], boxMax[d]);
    }
}


// The far-field model against the Python reference; the model's accuracy against a direct
// sum beyond 10 sigma
TEST(CoreFitSCTest, FarFieldAccuracy) {
    const Bunch b = makeSynthetic(false);
    const FarModel m = computeFarModel(b, syntheticGamma);

    // reference: exact core field plus the same far-far sum, so only the model is tested
    std::vector<Vector_t> reference(m.farFar);
    addDirectField(m.position, 0, m.position.size(), m.coreSource, m.coreCharge, 0.0, reference);
    double maxError = 0.0;
    std::size_t numBeyond = 0;
    for (std::size_t k = 0; k < m.table.size(); ++k) {
        const Vector_t dr = b.R[m.table[k]] - m.selection.mean;
        double dist = 0.0;
        for (unsigned int d = 0; d < 3; ++d) {
            dist = std::max(dist, std::abs(dr[d]) / m.selection.sigma[d]);
        }
        if (dist <= haloSigma) {
            continue;
        }
        ++ numBeyond;
        EXPECT_LE(m.indicator[k], indicatorThreshold);
        const Vector_t E = m.multipole[k] + m.farFar[k];
        maxError = std::max(maxError, euclidean_norm(E - reference[k]) /
                                      euclidean_norm(reference[k]));
    }
    EXPECT_EQ(numBeyond, 50u);
    EXPECT_LT(maxError, 5e-3);      // 1.8e-3 observed
}


// The model's accuracy close to the core. The 1 % tails of the scraped bunch put far
// particles just outside the range of cell centres, where the expansion is weakest. Those
// within 2 cells of the core bounds get the model (the exact core field where the indicator
// flags them), which is compared with the exact core field, both with the same far-far
// sum. The limits are floors against gross errors.
TEST(CoreFitSCTest, FarFieldAccuracyNearCore) {
    const Bunch b = makeSynthetic(true);
    const FarModel m = computeFarModel(b, syntheticGamma);

    std::vector<std::size_t> rows;
    std::vector<Vector_t> target;
    for (std::size_t k = 0; k < m.table.size(); ++k) {
        const Vector_t& r = b.R[m.table[k]];
        if (m.classes[m.table[k]] != ParticleClass::FAR) {
            continue;
        }
        double cells = 0.0;
        for (unsigned int d = 0; d < 3; ++d) {
            const double outside = std::max(m.selection.boundsMin[d] - r[d],
                                            r[d] - m.selection.boundsMax[d]);
            cells = std::max(cells, outside / m.hr[d]);
        }
        if (cells <= 2.0) {
            rows.push_back(k);
            target.push_back(m.position[k]);
        }
    }
    std::vector<Vector_t> exact(target.size(), Vector_t(0.0));
    addDirectField(target, 0, target.size(), m.coreSource, m.coreCharge, 0.0, exact);

    std::vector<double> error;
    std::size_t numFlagged = 0;
    for (std::size_t l = 0; l < rows.size(); ++l) {
        const std::size_t k = rows[l];
        const Vector_t reference = exact[l] + m.farFar[k];
        if (needsExactCoreField(m.classes[m.table[k]], m.indicator[k])) {
            ++ numFlagged;
            error.push_back(0.0);
        } else {
            const Vector_t E = m.multipole[k] + m.farFar[k];
            error.push_back(euclidean_norm(E - reference) / euclidean_norm(reference));
        }
    }
    ASSERT_GT(error.size(), 100u);  // 131, of them 21 flagged
    EXPECT_GT(numFlagged, 0u);
    std::sort(error.begin(), error.end());
    EXPECT_LT(error[error.size() / 2], 0.02);  // 2.5e-3 observed
    EXPECT_LT(error.back(), 0.08);              // 1.8e-2 observed
}


// The quadrupole term. Three point pairs +-z_k u about a centroid displaced from the box
// centre, with charges (1, -0.1, 1/135) q at z = (1, 2, 3) L: the l = 4 and l = 6 moments
// vanish, so beyond monopole and quadrupole the exact field holds only l >= 8 terms, about
// (L/r)^6 = 1e-18 of the quadrupole term at r = 1000 L. The quadrupole term alone (charge
// set to 0) must equal the exact field minus the monopole (long double) to 1e-12. A single
// pair leaves the l = 4 term, which must fall as (a/r)^2.
TEST(CoreFitSCTest, QuadrupoleTermOfDisplacedPointPairs) {
    const double L = 1e-3;
    const double q = 1e-15;
    const Vector_t u = Vector_t({1.0, 2.0, -0.5}) / euclidean_norm(Vector_t({1.0, 2.0, -0.5}));
    const Vector_t c({3.0 * L, -2.0 * L, 1.0 * L});

    std::vector<Vector_t> source;
    std::vector<double> charge;
    const double z[3] = {1.0, 2.0, 3.0};
    const double qk[3] = {q, -0.1 * q, q / 135.0};
    for (unsigned int k = 0; k < 3; ++k) {
        for (double sign: {1.0, -1.0}) {
            source.push_back(c + (sign * z[k] * L) * u);
            charge.push_back(qk[k]);
        }
    }
    MomentSums sums(Vector_t(0.0), 1.0);
    for (std::size_t j = 0; j < source.size(); ++j) {
        sums.add(charge[j], source[j]);
    }
    const CoreMoments moments = sums.reduce(SerialReducer());

    // moments: Q_ij = sum_k 2 q_k (z_k L)^2 (3 u_i u_j - delta_ij), centroid c
    double sz2 = 0.0;
    for (unsigned int k = 0; k < 3; ++k) {
        sz2 += 2.0 * qk[k] * z[k] * z[k] * L * L;
    }
    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            const double expected = sz2 * (3.0 * u[i] * u[j] - (i == j ? 1.0 : 0.0));
            EXPECT_LE(std::abs(moments.quadrupole[i][j] - expected), 1e-12 * 2.0 * sz2)
                << "Q_" << i << j;
        }
    }
    expectVectorNear(moments.centroid, c, 1e-12, "centroid");
    expectRelativeNear(moments.charge, 2.0 * (qk[0] + qk[1] + qk[2]), 1e-15, "charge");

    // the centroid in long double from the same positions, for the exact monopole
    long double qTotal = 0.0L;
    Vektor<long double, 3> centroid(0.0L);
    for (std::size_t j = 0; j < source.size(); ++j) {
        qTotal += charge[j];
        for (unsigned int k = 0; k < 3; ++k) {
            centroid[k] += charge[j] * static_cast<long double>(source[j][k]);
        }
    }
    for (unsigned int k = 0; k < 3; ++k) {
        centroid[k] /= qTotal;
    }

    CoreMoments quadrupoleOnly = moments;
    quadrupoleOnly.charge = 0.0;
    std::vector<Vector_t> target(1), field(1);
    std::vector<double> indicator(1);
    for (const Vector_t& dir: directions()) {
        target[0] = c + (1000.0 * L) * dir;
        computeMultipoleField(quadrupoleOnly, target, 0, 1, true, field, indicator);
        const Vektor<long double, 3> exact = exactField(target[0], source, charge);
        Vektor<long double, 3> d;
        for (unsigned int k = 0; k < 3; ++k) {
            d[k] = static_cast<long double>(target[0][k]) - centroid[k];
        }
        const long double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        const long double w = qTotal / (r2 * std::sqrt(r2));
        Vector_t reference;
        for (unsigned int k = 0; k < 3; ++k) {
            reference[k] = double(exact[k] - w * d[k]);
        }
        expectVectorNear(field[0], reference, 1e-12, "quadrupole term");
    }

    // a single displaced pair: the residual after monopole + quadrupole is the l = 4 term
    std::vector<Vector_t> pair = {c + L * u, c - L * u};
    std::vector<double> pairCharge = {q, q};
    MomentSums pairSums(Vector_t(0.0), 1.0);
    pairSums.add(q, pair[0]);
    pairSums.add(q, pair[1]);
    const CoreMoments pairMoments = pairSums.reduce(SerialReducer());
    const Vector_t dir = Vector_t({0.3, -0.7, 0.2}) / euclidean_norm(Vector_t({0.3, -0.7, 0.2}));
    double residual[2];
    for (unsigned int k = 0; k < 2; ++k) {
        target[0] = c + ((100.0 * (k + 1)) * L) * dir;
        computeMultipoleField(pairMoments, target, 0, 1, true, field, indicator);
        const Vektor<long double, 3> exact = exactField(target[0], pair, pairCharge);
        Vector_t diff;
        for (unsigned int d = 0; d < 3; ++d) {
            diff[d] = double(exact[d] - static_cast<long double>(field[0][d]));
        }
        CoreMoments quadOnly = pairMoments;
        quadOnly.charge = 0.0;
        std::vector<Vector_t> quadField(1);
        computeMultipoleField(quadOnly, target, 0, 1, true, quadField, indicator);
        residual[k] = euclidean_norm(diff) / euclidean_norm(quadField[0]);
    }
    EXPECT_LT(residual[0], 1e-3);
    EXPECT_NEAR(residual[0] / residual[1], 4.0, 0.05);
}


// The rest frame. The model with rest-frame positions and the lab transform reproduces the
// exact field of a point charge in uniform motion along y,
// E = Q (1 - beta^2) R / (|R|^3 (1 - beta^2 sin^2 theta)^(3/2)); with lab positions the
// longitudinal field is wrong by gamma^2 - 1.
TEST(CoreFitSCTest, RestFrameAndLabTransform) {
    const double Q = 1e-15;
    const Vector_t charge({1.2e-3, -0.4e-3, 0.3e-3});   // local frame
    const Vector_t centre({1.0e-3, 0.5e-3, -0.2e-3});   // box centre
    for (double gamma: {syntheticGamma, 3.0}) {
        const double beta2 = 1.0 - 1.0 / (gamma * gamma);
        MomentSums sums(centre, gamma);
        sums.add(Q, charge);
        const CoreMoments moments = sums.reduce(SerialReducer());
        expectVectorNear(moments.centroid, toRestFrame(charge, centre, gamma), 1e-15, "centroid");

        MomentSums labSums(centre, 1.0);
        labSums.add(Q, charge);
        const CoreMoments labMoments = labSums.reduce(SerialReducer());

        std::vector<Vector_t> target(1), field(1);
        std::vector<double> indicator(1);
        for (const Vector_t& dir: directions()) {
            const Vector_t R = 0.05 * dir;
            const Vector_t r = charge + R;
            const double sin2 = (R[0] * R[0] + R[2] * R[2]) / dot(R, R);
            const double Rn = euclidean_norm(R);
            const Vector_t exact = (Q * (1.0 - beta2) /
                                    (Rn * Rn * Rn * std::pow(1.0 - beta2 * sin2, 1.5))) * R;

            target[0] = toRestFrame(r, centre, gamma);
            computeMultipoleField(moments, target, 0, 1, true, field, indicator);
            expectVectorNear(toLabFrame(field[0], gamma), exact, 1e-13, "rest frame");
            EXPECT_LE(indicator[0], 1e-10);

            // lab positions: the pitfall
            target[0] = r - centre;
            computeMultipoleField(labMoments, target, 0, 1, true, field, indicator);
            const Vector_t wrong = toLabFrame(field[0], gamma);
            if (std::abs(dir[1]) == 1.0) {
                expectRelativeNear(wrong[1] / exact[1] - 1.0, gamma * gamma - 1.0, 1e-12,
                                   "lab positions, longitudinal");
            }
            if (dir[1] == 0.0) {
                expectVectorNear(wrong, exact, 1e-13, "lab positions, transverse");
            }
        }
    }
}


// The work split over P = 1, 3 and 8 ranks gives bitwise the same fields, and each rank
// gets n / P rows of the table, within one. The lost entries would add nothing as far-far
// sources, and a pair at the same position is skipped.
TEST(CoreFitSCTest, WorkSplitByTableIndex) {
    const Bunch b = makeSynthetic(true);
    const FarModel m = computeFarModel(b, syntheticGamma);
    const std::size_t n = m.table.size();
    ASSERT_GT(n, 500u);

    for (int numRanks: {1, 3, 8}) {
        std::vector<Vector_t> multipole(n, Vector_t(0.0)), farFar(n, Vector_t(0.0));
        std::vector<double> indicator(n, -1.0);
        std::size_t covered = 0;
        for (int p = 0; p < numRanks; ++p) {
            const std::pair<std::size_t, std::size_t> range = getTableRange(n, p, numRanks);
            EXPECT_EQ(range.first, covered);
            covered = range.second;
            computeMultipoleField(m.moments, m.position, range.first, range.second, true,
                                  multipole, indicator);
            addDirectField(m.position, range.first, range.second, m.farSource, m.farCharge,
                           m.eps, farFar);
            const double rows = double(range.second - range.first);
            EXPECT_LE(std::abs(rows - double(n) / numRanks), 1.0);
        }
        EXPECT_EQ(covered, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (unsigned int d = 0; d < 3; ++d) {
                EXPECT_EQ(multipole[i][d], m.multipole[i][d]);
                EXPECT_EQ(farFar[i][d], m.farFar[i][d]);
            }
            EXPECT_EQ(indicator[i], m.indicator[i]);
        }
    }
    EXPECT_EQ(getTableRange(10, 1, 3), std::make_pair(std::size_t(3), std::size_t(6)));
    // more ranks than entries
    std::size_t covered = 0;
    for (int p = 0; p < 8; ++p) {
        const std::pair<std::size_t, std::size_t> range = getTableRange(3, p, 8);
        EXPECT_EQ(range.first, covered);
        EXPECT_LE(range.second - range.first, 1u);
        covered = range.second;
    }
    EXPECT_EQ(covered, 3u);

    // the monopole alone
    std::vector<Vector_t> monopole(n);
    std::vector<double> indicator(n);
    computeMultipoleField(m.moments, m.position, 0, n, false, monopole, indicator);
    for (std::size_t i = 0; i < n; ++i) {
        expectRelativeNear(indicator[i], m.indicator[i], 1e-14, "indicator");
        const Vector_t d = m.position[i] - m.moments.centroid;
        expectVectorNear(monopole[i], (m.moments.charge / std::pow(dot(d, d), 1.5)) * d,
                         1e-15, "monopole");
    }

    // the lost entries (charge 0) as extra sources change no bit of the far-far field
    std::vector<Vector_t> withLost(n, Vector_t(0.0));
    addDirectField(m.position, 0, n, m.position, m.charge, m.eps, withLost);
    ASSERT_GT(m.position.size(), m.farSource.size());
    for (std::size_t i = 0; i < n; ++i) {
        for (unsigned int d = 0; d < 3; ++d) {
            EXPECT_EQ(withLost[i][d], m.farFar[i][d]);
        }
    }

    // a target on a source: that pair is skipped, also without softening
    const std::vector<Vector_t> points = {Vector_t({1e-3, 0.0, 0.0}), Vector_t({0.0, 2e-3, 0.0}),
                                          Vector_t({0.0, 0.0, -1e-3})};
    const std::vector<double> charges = {1e-15, 2e-15, 3e-15};
    std::vector<Vector_t> onSource(3, Vector_t(0.0));
    addDirectField(points, 0, 3, points, charges, 0.0, onSource);
    for (std::size_t k = 0; k < 3; ++k) {
        std::vector<Vector_t> others;
        std::vector<double> otherCharges;
        for (std::size_t j = 0; j < 3; ++j) {
            if (j != k) {
                others.push_back(points[j]);
                otherCharges.push_back(charges[j]);
            }
        }
        const Vektor<long double, 3> exact = exactField(points[k], others, otherCharges);
        expectVectorNear(onSource[k], Vector_t({double(exact[0]), double(exact[1]),
                                                double(exact[2])}), 1e-15, "target on a source");
    }
}


// The subsample of the sources of a direct sum above its pair budget. The stride keeps the
// pairs within the budget. Every ID is in the subsample of stride 1, and about one in
// stride otherwise, for consecutive IDs and for IDs in steps of the stride (IPPL numbers the
// particles created on one rank in steps of the rank count). The outer shell of the box is
// the part of it beyond 1 - coreShellFraction of its half-width in some axis. On the scraped
// bunch, the exact core field at the flagged FAR entries with the budget of the solve, the
// core particles in the shell and the others each subsampled to maxCorePairsPerCore pairs
// per core particle and scaled to their charge, is close to the field of all core
// particles, and much closer than the multipole field that the indicator rejects there.
TEST(CoreFitSCTest, DirectSumSubsample) {
    EXPECT_EQ(getSampleStride(0, 0, maxFarFarPairs), 1u);
    EXPECT_EQ(getSampleStride(10000, 10000, maxFarFarPairs), 1u);
    EXPECT_EQ(getSampleStride(10000, 10001, maxFarFarPairs), 2u);
    EXPECT_EQ(getSampleStride(20000, 20000, maxFarFarPairs), 4u);
    // the exact core field: about one in N_flag / maxCorePairsPerCore core particles
    EXPECT_EQ(getSampleStride(8, 349641, maxCorePairsPerCore * 349641), 1u);
    EXPECT_EQ(getSampleStride(9, 349641, maxCorePairsPerCore * 349641), 2u);
    EXPECT_EQ(getSampleStride(292, 349641, maxCorePairsPerCore * 349641), 37u);

    const Vector_t lower({-1.0, -2.0, -3.0});
    const Vector_t upper({1.0, 2.0, 5.0});
    EXPECT_FALSE(isInCoreShell(Vector_t({0.0, 0.0, 1.0}), lower, upper));
    EXPECT_FALSE(isInCoreShell(Vector_t({0.74, -1.49, 3.99}), lower, upper));
    EXPECT_TRUE(isInCoreShell(Vector_t({0.76, 0.0, 1.0}), lower, upper));
    EXPECT_TRUE(isInCoreShell(Vector_t({0.0, -1.51, 1.0}), lower, upper));
    EXPECT_TRUE(isInCoreShell(Vector_t({0.0, 0.0, -2.01}), lower, upper));
    EXPECT_TRUE(isInCoreShell(lower, lower, upper));
    EXPECT_TRUE(isInCoreShell(upper, lower, upper));

    for (std::size_t id = 0; id < 1000; ++id) {
        EXPECT_TRUE(isSampled(id, 1));
    }
    const std::size_t numIds = 100000;
    for (std::size_t stride: {2, 3, 8, 100}) {
        std::size_t consecutive = 0;
        std::size_t steps = 0;
        for (std::size_t id = 0; id < numIds; ++id) {
            consecutive += isSampled(id, stride);
            steps += isSampled(id * stride, stride);
        }
        const double expected = double(numIds) / stride;
        EXPECT_LT(std::abs(consecutive - expected), 5.0 * std::sqrt(expected));
        EXPECT_LT(std::abs(steps - expected), 5.0 * std::sqrt(expected));
    }

    const Bunch b = makeSynthetic(true);
    const FarModel m = computeFarModel(b, syntheticGamma);
    ASSERT_GT(m.flagged.size(), 50u);
    std::vector<Vector_t> target;
    for (std::size_t k: m.flagged) {
        target.push_back(m.position[k]);
    }
    std::vector<Vector_t> exact(target.size(), Vector_t(0.0));
    addDirectField(target, 0, target.size(), m.coreSource, m.coreCharge, 0.0, exact);

    // the core particles in index order, with their index as ID, in the shell of the box
    // (the core bounds) or not
    const std::size_t numCore = m.coreSource.size();
    const double maxPairs = maxCorePairsPerCore * numCore;
    std::vector<bool> inShell;
    std::size_t numShell = 0;
    for (std::size_t i = 0; i < b.R.size(); ++i) {
        if (m.classes[i] == ParticleClass::CORE) {
            inShell.push_back(isInCoreShell(b.R[i], m.selection.boundsMin,
                                            m.selection.boundsMax));
            numShell += inShell.back();
        }
    }
    const std::size_t stride[2] = {getSampleStride(target.size(), numShell, maxPairs),
                                   getSampleStride(target.size(), numCore - numShell,
                                                   maxPairs)};
    // the numpy implementation: 475 of the 201241 core particles in the shell, strides 1 and 12
    EXPECT_GT(numShell, 100u);
    EXPECT_LT(numShell, numCore / 100);
    EXPECT_EQ(stride[0], 1u);
    EXPECT_GT(stride[1], 8u);

    std::vector<Vector_t> sampled(target.size(), Vector_t(0.0));
    std::size_t numSources = 0;
    for (int part = 0; part < 2; ++part) {
        std::vector<Vector_t> source;
        std::vector<double> charge;
        double partCharge = 0.0;
        double sampledCharge = 0.0;
        for (std::size_t i = 0, l = 0; i < b.R.size(); ++i) {
            if (m.classes[i] != ParticleClass::CORE) {
                continue;
            }
            if (inShell[l] == (part == 0)) {
                partCharge += m.coreCharge[l];
                if (isSampled(i, stride[part])) {
                    source.push_back(m.coreSource[l]);
                    charge.push_back(m.coreCharge[l]);
                    sampledCharge += m.coreCharge[l];
                }
            }
            ++ l;
        }
        for (double& q: charge) {
            q *= partCharge / sampledCharge;
        }
        addDirectField(target, 0, target.size(), source, charge, 0.0, sampled);
        numSources += source.size();
    }
    EXPECT_LT(std::abs(double(numSources) - double(numShell) -
                       double(numCore - numShell) / stride[1]),
              5.0 * std::sqrt(double(numCore - numShell) / stride[1]));
    EXPECT_LE(double(numSources) * target.size(), 2.0 * maxPairs);

    std::vector<double> error, multipoleError;
    for (std::size_t l = 0; l < target.size(); ++l) {
        const double norm = euclidean_norm(exact[l]);
        error.push_back(euclidean_norm(sampled[l] - exact[l]) / norm);
        multipoleError.push_back(euclidean_norm(m.multipole[m.flagged[l]] - exact[l]) / norm);
    }
    std::sort(error.begin(), error.end());
    std::sort(multipoleError.begin(), multipoleError.end());
    const std::size_t median = error.size() / 2;
    // the numpy implementation gives 4.9e-3 and 9.0e-3 on the same particles, the multipole
    // 6.0e-2 and 0.39
    EXPECT_LT(error[median], 1e-2);
    EXPECT_LT(error.back(), 2e-2);
    EXPECT_LT(error[median], 0.2 * multipoleError[median]);
    EXPECT_LT(error.back(), 0.2 * multipoleError.back());
}


TEST(CoreFitSCTest, Classification) {
    const Vector_t boundsMin(-1.0), boundsMax(1.0), meshMin(-1.1), meshMax(1.1);
    auto cls = [&](int bin, const Vector_t& r) {
        return classify(bin, r, boundsMin, boundsMax, meshMin, meshMax);
    };
    EXPECT_EQ(cls(0, Vector_t(0.0)), ParticleClass::CORE);
    EXPECT_EQ(cls(3, Vector_t({1.0, -1.0, 1.0})), ParticleClass::CORE);
    EXPECT_EQ(cls(-1, Vector_t(0.5)), ParticleClass::INSIDE);
    EXPECT_EQ(cls(0, Vector_t({1.05, 0.0, 0.0})), ParticleClass::MESH);
    EXPECT_EQ(cls(-1, Vector_t({0.0, -1.1, 0.0})), ParticleClass::MESH);
    EXPECT_EQ(cls(0, Vector_t({0.0, 0.0, 1.1000001})), ParticleClass::FAR);
    EXPECT_EQ(cls(-1, Vector_t({-5.0, 0.0, 0.0})), ParticleClass::FAR);
    EXPECT_EQ(cls(0, Vector_t({1.05, 1.2, 0.0})), ParticleClass::FAR);

    // the exact core field replaces the multipole field of a FAR entry above the threshold
    // only; the MESH band is gathered from the mesh whatever its indicator
    EXPECT_TRUE(needsExactCoreField(ParticleClass::FAR, 0.2));
    EXPECT_TRUE(needsExactCoreField(ParticleClass::FAR,
                                    std::nextafter(indicatorThreshold, 1.0)));
    EXPECT_FALSE(needsExactCoreField(ParticleClass::FAR, indicatorThreshold));
    EXPECT_FALSE(needsExactCoreField(ParticleClass::FAR, 0.0));
    for (ParticleClass c: {ParticleClass::CORE, ParticleClass::INSIDE, ParticleClass::MESH}) {
        EXPECT_FALSE(needsExactCoreField(c, 0.2));
        EXPECT_FALSE(needsExactCoreField(c, 1e3));
    }
}


// The mode decision with hysteresis
TEST(CoreFitSCTest, ModeDecision) {
    const unsigned int minDwell = 10;       // REPARTFREQ
    const double maxFar = 0.05;             // MESHFITMAXFAR
    const std::size_t live = 200000;
    const std::size_t core = 199000;

    // a run starts from ModeState{} = {CORE, 0}: the first solve may switch, and a solve
    // without a switch counts as the minimum dwell reached
    const ModeState start;
    EXPECT_EQ(start.mode, Mode::CORE);
    EXPECT_EQ(start.dwell, 0u);
    ModeState next = decideMode(start, core, live, 0.06, maxFar, minDwell);
    EXPECT_EQ(next.mode, Mode::FULL);
    EXPECT_EQ(next.dwell, 1u);
    next = decideMode(start, core, live, 0.0, maxFar, minDwell);
    EXPECT_EQ(next.mode, Mode::CORE);
    EXPECT_EQ(next.dwell, minDwell);
    EXPECT_EQ(decideMode(next, core, live, 0.06, maxFar, minDwell).mode, Mode::FULL);
    // the same from {CORE, minDwell}
    ModeState state{Mode::CORE, minDwell};
    next = decideMode(state, core, live, 0.06, maxFar, minDwell);
    EXPECT_EQ(next.mode, Mode::FULL);
    EXPECT_EQ(next.dwell, 1u);

    // a far fraction oscillating around the limit: one switch to FULL, and no return while
    // it stays above half of the limit
    state = ModeState{Mode::CORE, minDwell};
    unsigned int switches = 0;
    for (int step = 0; step < 200; ++step) {
        const double far = (step % 2 == 0 ? 0.04 : 0.06);
        next = decideMode(state, core, live, far, maxFar, minDwell);
        if (next.mode != state.mode) {
            ++ switches;
            EXPECT_EQ(step, 1);
        }
        state = next;
    }
    EXPECT_EQ(switches, 1u);
    EXPECT_EQ(state.mode, Mode::FULL);
    EXPECT_EQ(state.dwell, minDwell);

    // oscillating between 0.02 and 0.06: at most one switch per minDwell steps, and CORE is
    // entered only on a step below half of the limit
    state = ModeState{Mode::CORE, minDwell};
    int lastSwitch = -1000;
    switches = 0;
    for (int step = 0; step < 300; ++step) {
        const double far = (step % 3 == 0 ? 0.02 : 0.06);
        next = decideMode(state, core, live, far, maxFar, minDwell);
        if (next.mode != state.mode) {
            ++ switches;
            EXPECT_GE(step - lastSwitch, int(minDwell));
            if (next.mode == Mode::CORE) {
                EXPECT_LT(far, 0.5 * maxFar);
            } else {
                EXPECT_GT(far, maxFar);
            }
            lastSwitch = step;
        }
        state = next;
    }
    EXPECT_GT(switches, 10u);

    // no switch before the dwell is reached
    state = ModeState{Mode::FULL, 1};
    for (unsigned int step = 1; step < minDwell; ++step) {
        state = decideMode(state, core, live, 0.0, maxFar, minDwell);
        EXPECT_EQ(state.mode, Mode::FULL);
        EXPECT_EQ(state.dwell, step + 1);
    }
    state = decideMode(state, core, live, 0.0, maxFar, minDwell);
    EXPECT_EQ(state.mode, Mode::CORE);

    // ... except for a broken-up or empty core, which leaves CORE at once; a large far
    // fraction alone waits for the dwell
    for (unsigned int dwell = 1; dwell < minDwell; ++dwell) {
        state = ModeState{Mode::CORE, dwell};
        EXPECT_EQ(decideMode(state, 0, live, 1.0, maxFar, minDwell).mode, Mode::FULL);
        EXPECT_EQ(decideMode(state, 500, live, 0.9, maxFar, minDwell).mode, Mode::FULL);
        EXPECT_EQ(decideMode(state, 99999, live, 0.0, maxFar, minDwell).mode, Mode::FULL);
        next = decideMode(state, core, live, 0.06, maxFar, minDwell);
        EXPECT_EQ(next.mode, Mode::CORE);
        EXPECT_EQ(next.dwell, dwell + 1);
    }

    // the core count: FULL below 1000 or below half of the live particles
    state = ModeState{Mode::CORE, minDwell};
    EXPECT_EQ(decideMode(state, 999, 1500, 0.0, maxFar, minDwell).mode, Mode::FULL);
    EXPECT_EQ(decideMode(state, 1000, 1500, 0.0, maxFar, minDwell).mode, Mode::CORE);
    EXPECT_EQ(decideMode(state, 99999, 200000, 0.0, maxFar, minDwell).mode, Mode::FULL);
    EXPECT_EQ(decideMode(state, 100000, 200000, 0.0, maxFar, minDwell).mode, Mode::CORE);
    EXPECT_EQ(decideMode(state, 0, 0, 0.0, maxFar, minDwell).mode, Mode::FULL);
    // back to CORE only above 60 % and at least 1000
    state = ModeState{Mode::FULL, minDwell};
    EXPECT_EQ(decideMode(state, 120000, 200000, 0.0, maxFar, minDwell).mode, Mode::FULL);
    EXPECT_EQ(decideMode(state, 120001, 200000, 0.0, maxFar, minDwell).mode, Mode::CORE);
    EXPECT_EQ(decideMode(state, 999, 1200, 0.0, maxFar, minDwell).mode, Mode::FULL);
    EXPECT_EQ(decideMode(state, 1000, 1200, 0.0, maxFar, minDwell).mode, Mode::CORE);
    EXPECT_EQ(decideMode(state, 199000, 200000, 0.025, maxFar, minDwell).mode, Mode::FULL);

    // REPARTFREQ = 1 (and 0): a switch is allowed on every step, and the dwell saturates at 1
    EXPECT_EQ(decideMode(ModeState{}, core, live, 0.0, maxFar, 0).dwell, 1u);
    for (unsigned int dwell: {0u, 1u}) {
        state = ModeState{Mode::CORE, dwell};
        for (int step = 0; step < 6; ++step) {
            const double far = (step % 2 == 0 ? 0.06 : 0.0);
            next = decideMode(state, core, live, far, maxFar, dwell);
            EXPECT_NE(next.mode, state.mode) << "step " << step;
            state = next;
        }
    }
}
