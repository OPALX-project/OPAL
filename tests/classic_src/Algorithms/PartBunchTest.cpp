//
// Tests for PartBunch::scatterMasked() and PartBunch::gatherMasked()
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
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//

#include "gtest/gtest.h"

#include "Algorithms/PartBunch.h"
#include "Algorithms/PartData.h"
#include "Physics/Physics.h"

#include "opal_test_utilities/SilenceTest.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <random>
#include <vector>

namespace {
    // Cells per axis. With four ranks every vnode is three cells wide along z.
    constexpr int nCells = 12;
    // Particles inside the mesh, all created on rank 0 and distributed by update().
    constexpr size_t nInside = 4000;
    // Particles outside the mesh, used on one rank only: without a layout policy for
    // them, update() cannot place them on more than one rank.
    constexpr size_t nOutside = 4;

    class PartBunchTest: public ::testing::Test {
    public:
        PartBunchTest():
            mesh_m(meshDomain()),
            fieldLayout_m(mesh_m, decomp_m),
            partData_m(1.0, Physics::m_p * 1e9, 1e8),
            bunch_m(&partData_m),
            rmin_m({-0.31, 1.07, -0.012}),
            rmax_m({-0.26, 1.13, 0.009}),
            sentinel_m({1.5e30, -2.5e30, 3.5e30})
        { }

        void SetUp() override {
            bunch_m.initialize(&fieldLayout_m);

            // Mesh as in PartBunchBase::boundp(): nCells cell centres from rmin to rmax.
            for (unsigned d = 0; d < 3; ++d) {
                hr_m(d) = (rmax_m(d) - rmin_m(d)) / (nCells - 1);
            }
            Vector_t origin = rmin_m - 0.5 * hr_m;
            bunch_m.getMesh().set_meshSpacing(&(hr_m[0]));
            bunch_m.getMesh().set_origin(origin);

            for (int i = 0; i < 2 * 3; ++i) {
                bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
                vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
            }
            bunch_m.rho_m.initialize(bunch_m.getMesh(), bunch_m.getFieldLayout(),
                                     GuardCellSizes<3>(1), bc_m);
            bunch_m.eg_m.initialize(bunch_m.getMesh(), bunch_m.getFieldLayout(),
                                    GuardCellSizes<3>(1), vbc_m);
        }

    protected:
        static NDIndex<3> meshDomain() {
            NDIndex<3> domain;
            for (unsigned d = 0; d < 3; ++d) {
                domain[d] = Index(nCells + 1);
            }
            return domain;
        }

        // Uniform in [rmin, rmax], plus two particles on the corner cell centres; every
        // seventh particle is flagged lost. With withOutside, nOutside more particles lie
        // outside the mesh: less than a cell outside (in the guard layer) and far outside.
        void createParticles(bool withOutside) {
            const size_t nOut = withOutside ? nOutside : 0;
            if (Ippl::myNode() == 0) {
                std::mt19937_64 generator(20260924);
                std::uniform_real_distribution<double> unit(0.0, 1.0);
                bunch_m.create(nInside + nOut);
                for (size_t i = 0; i < nInside + nOut; ++i) {
                    for (unsigned d = 0; d < 3; ++d) {
                        bunch_m.R[i](d) = rmin_m(d) + unit(generator) * (rmax_m(d) - rmin_m(d));
                    }
                    bunch_m.Q[i] = 1e-15 * (0.5 + unit(generator));
                    bunch_m.Bin[i] = (i % 7 == 3) ? -1 : 0;
                }
                bunch_m.R[0] = rmin_m;
                bunch_m.R[1] = rmax_m;
                for (size_t i = nInside; i < nInside + nOut; ++i) {
                    const size_t k = i - nInside;
                    const unsigned d = k % 3;
                    const double shift = (k % 2 == 0) ? 1.2 : 40.0;
                    bunch_m.R[i](d) = (k < 2) ? rmin_m(d) - shift * hr_m(d)
                                              : rmax_m(d) + shift * hr_m(d);
                }
            }
            bunch_m.update();
            ASSERT_EQ(bunch_m.getTotalNum(), nInside + nOut);
            // every rank owns a slab of the mesh and some of the particles
            EXPECT_GT(bunch_m.getLocalNum(), 0u);
            if (Ippl::getNodes() > 1) {
                EXPECT_LT(bunch_m.getLocalNum(), nInside);
            }
            // and the particles inside the mesh are those of its own slab, which the CIC
            // kernels need
            const NDIndex<3> local = bunch_m.getFieldLayout().getLocalNDIndex();
            const double zlo = bunch_m.getMesh().get_origin()(2) + local[2].first() * hr_m(2);
            const double zhi = zlo + local[2].length() * hr_m(2);
            size_t nForeign = 0;
            for (size_t i = 0; i < bunch_m.getLocalNum(); ++i) {
                const double z = bunch_m.R[i](2);
                if (insideCentres(bunch_m.R[i]) && (z < zlo || z > zhi)) {
                    ++nForeign;
                }
            }
            EXPECT_EQ(globalSum(nForeign), 0u);
        }

        bool insideCentres(const Vector_t& x) const {
            for (unsigned d = 0; d < 3; ++d) {
                if (x(d) < rmin_m(d) || x(d) > rmax_m(d)) {
                    return false;
                }
            }
            return true;
        }

        // A predicate that rejects about a third of the particles inside the mesh and
        // every particle outside it.
        PartBunch::ParticleMask_t partialMask() const {
            return [this](size_t i) {
                return insideCentres(bunch_m.R[i]) && bunch_m.ID[i] % 3 != 0;
            };
        }

        // The scatter of computeSelfFields_cycl(double), with the charge of the particles
        // rejected by accept set to zero and those outside the mesh moved onto a cell centre
        // for the call.
        std::vector<double> referenceScatter(const PartBunch::ParticleMask_t& accept) {
            const size_t localNum = bunch_m.getLocalNum();
            std::vector<double> q(localNum);
            std::vector<Vector_t> r(localNum);
            for (size_t i = 0; i < localNum; ++i) {
                q[i] = bunch_m.Q[i];
                r[i] = bunch_m.R[i];
                if (bunch_m.Bin[i] < 0 || (accept && !accept(i))) {
                    bunch_m.Q[i] = 0.0;
                }
                if (!insideCentres(bunch_m.R[i])) {
                    bunch_m.R[i] = rmin_m;
                }
            }
            bunch_m.rho_m = 0.0;
            bunch_m.Q.scatter(bunch_m.rho_m, bunch_m.R, IntrplCIC_t());
            for (size_t i = 0; i < localNum; ++i) {
                bunch_m.Q[i] = q[i];
                bunch_m.R[i] = r[i];
            }
            return localValues(bunch_m.rho_m);
        }

        // Ef.gather() on every particle, with those outside the mesh moved onto a cell
        // centre for the call (their values are not used).
        std::vector<Vector_t> referenceGather() {
            const size_t localNum = bunch_m.getLocalNum();
            std::vector<Vector_t> r(localNum);
            for (size_t i = 0; i < localNum; ++i) {
                r[i] = bunch_m.R[i];
                if (!insideCentres(bunch_m.R[i])) {
                    bunch_m.R[i] = rmin_m;
                }
            }
            setEf(sentinel_m);
            bunch_m.Ef.gather(bunch_m.eg_m, bunch_m.R, IntrplCIC_t());
            for (size_t i = 0; i < localNum; ++i) {
                bunch_m.R[i] = r[i];
            }
            return localEf();
        }

        // Some field on eg_m: the gradient of the charge of the particles inside the mesh.
        // With deferred guard cell fills (--defergcfill) eg_m is left dirty: its internal
        // guard cells hold values from before the last steps of Grad() and the negation,
        // which a gather must not use. They are set to NaN first, so that no value of an
        // earlier field is left in them either.
        void fillFieldFromCharge() {
            bunch_m.eg_m.setGuardCells(Vector_t(std::numeric_limits<double>::quiet_NaN()));
            bunch_m.rho_m = 0.0;
            bunch_m.scatterMasked([this](size_t i) { return insideCentres(bunch_m.R[i]); });
            bunch_m.rho_m *= 1e15;
            bunch_m.eg_m = -Grad(bunch_m.rho_m, bunch_m.eg_m);
        }

        void setEf(const Vector_t& value) {
            for (size_t i = 0; i < bunch_m.getLocalNum(); ++i) {
                bunch_m.Ef[i] = value;
            }
        }

        std::vector<Vector_t> localEf() const {
            std::vector<Vector_t> ef(bunch_m.getLocalNum());
            for (size_t i = 0; i < ef.size(); ++i) {
                ef[i] = bunch_m.Ef[i];
            }
            return ef;
        }

        std::vector<double> localQ() const {
            std::vector<double> q(bunch_m.getLocalNum());
            for (size_t i = 0; i < q.size(); ++i) {
                q[i] = bunch_m.Q[i];
            }
            return q;
        }

        // The owned cells of this rank, in iteration order.
        static std::vector<double> localValues(const Field_t& field) {
            std::vector<double> values;
            for (Field_t::iterator it = field.begin(); it != field.end(); ++it) {
                values.push_back(*it);
            }
            return values;
        }

        static bool bitIdentical(const double& a, const double& b) {
            return std::memcmp(&a, &b, sizeof(double)) == 0;
        }

        static bool bitIdentical(const Vector_t& a, const Vector_t& b) {
            for (unsigned d = 0; d < 3; ++d) {
                const double ad = a(d);
                const double bd = b(d);
                if (!bitIdentical(ad, bd)) {
                    return false;
                }
            }
            return true;
        }

        // Number of differing entries, summed over all ranks.
        template <class T>
        static size_t countDifferent(const std::vector<T>& a, const std::vector<T>& b) {
            size_t n = (a.size() == b.size()) ? 0 : 1;
            for (size_t i = 0; i < std::min(a.size(), b.size()); ++i) {
                if (!bitIdentical(a[i], b[i])) {
                    ++n;
                }
            }
            allreduce(n, 1, std::plus<size_t>());
            return n;
        }

        static size_t globalSum(size_t n) {
            allreduce(n, 1, std::plus<size_t>());
            return n;
        }

        // Number of entries with a component that is not finite, summed over all ranks.
        static size_t countNonFinite(const std::vector<Vector_t>& values) {
            size_t n = 0;
            for (const Vector_t& v : values) {
                if (!std::isfinite(v(0)) || !std::isfinite(v(1)) || !std::isfinite(v(2))) {
                    ++n;
                }
            }
            return globalSum(n);
        }

        e_dim_tag decomp_m[3] = {SERIAL, SERIAL, PARALLEL};
        Mesh_t mesh_m;
        FieldLayout_t fieldLayout_m;
        BConds<double, 3, Mesh_t, Center_t> bc_m;
        BConds<Vector_t, 3, Mesh_t, Center_t> vbc_m;
        PartData partData_m;
        PartBunch bunch_m;
        const Vector_t rmin_m;
        const Vector_t rmax_m;
        const Vector_t sentinel_m;
        Vector_t hr_m;
    };
}

TEST_F(PartBunchTest, ScatterMaskedAllTrueIsBitIdentical) {
    OpalTestUtilities::SilenceTest silencer;

    createParticles(false);
    const std::vector<double> expected = referenceScatter(PartBunch::ParticleMask_t());
    const std::vector<double> q = localQ();
    ASSERT_GT(globalSum(expected.size()), 0u);

    bunch_m.rho_m = 0.0;
    bunch_m.scatterMasked([](size_t) { return true; });
    EXPECT_EQ(countDifferent(localValues(bunch_m.rho_m), expected), 0u);
    EXPECT_EQ(countDifferent(localQ(), q), 0u);

    // and without a predicate
    bunch_m.rho_m = 0.0;
    bunch_m.scatterMasked(PartBunch::ParticleMask_t());
    EXPECT_EQ(countDifferent(localValues(bunch_m.rho_m), expected), 0u);
    EXPECT_EQ(countDifferent(localQ(), q), 0u);
}

TEST_F(PartBunchTest, ScatterMaskedLeavesOutRejected) {
    OpalTestUtilities::SilenceTest silencer;

    createParticles(Ippl::getNodes() == 1);
    const PartBunch::ParticleMask_t accept = partialMask();
    const std::vector<double> q = localQ();

    double acceptedCharge = 0.0;
    size_t nRejected = 0;
    for (size_t i = 0; i < bunch_m.getLocalNum(); ++i) {
        if (bunch_m.Bin[i] >= 0 && accept(i)) {
            acceptedCharge += bunch_m.Q[i];
        } else if (bunch_m.Bin[i] >= 0) {
            ++nRejected;
        }
    }
    allreduce(acceptedCharge, 1, std::plus<double>());
    EXPECT_GT(globalSum(nRejected), nInside / 5);

    bunch_m.rho_m = 0.0;
    bunch_m.scatterMasked(accept);
    const std::vector<double> masked = localValues(bunch_m.rho_m);
    const double meshCharge = sum(bunch_m.rho_m);
    EXPECT_NEAR(meshCharge, acceptedCharge, 1e-14 * acceptedCharge);
    EXPECT_EQ(countDifferent(localQ(), q), 0u);

    EXPECT_EQ(countDifferent(masked, referenceScatter(accept)), 0u);
    // and the comparison does see the rejected charge
    EXPECT_GT(countDifferent(masked, referenceScatter(PartBunch::ParticleMask_t())), 0u);

    // a predicate that accepts more particles than the one before
    const PartBunch::ParticleMask_t inside = [this](size_t i) {
        return insideCentres(bunch_m.R[i]);
    };
    bunch_m.rho_m = 0.0;
    bunch_m.scatterMasked(inside);
    EXPECT_EQ(countDifferent(localValues(bunch_m.rho_m), referenceScatter(inside)), 0u);
    EXPECT_EQ(countDifferent(localQ(), q), 0u);
}

TEST_F(PartBunchTest, GatherMaskedAllTrueIsBitIdentical) {
    OpalTestUtilities::SilenceTest silencer;

    createParticles(false);
    // eg_m is computed again before each gather, so that with deferred guard cell fills
    // every call starts from a dirty field
    fillFieldFromCharge();
    const std::vector<Vector_t> expected = referenceGather();
    size_t nNonzero = 0;
    for (const Vector_t& e : expected) {
        if (e(0) != 0.0 && e(1) != 0.0 && e(2) != 0.0) {
            ++nNonzero;
        }
    }
    EXPECT_GT(globalSum(nNonzero), nInside / 2);
    EXPECT_EQ(countNonFinite(expected), 0u);

    fillFieldFromCharge();
    setEf(sentinel_m);
    bunch_m.gatherMasked([](size_t) { return true; });
    EXPECT_EQ(countDifferent(localEf(), expected), 0u);
}

TEST_F(PartBunchTest, GatherMaskedLeavesRejectedUntouched) {
    OpalTestUtilities::SilenceTest silencer;

    createParticles(Ippl::getNodes() == 1);
    const PartBunch::ParticleMask_t accept = partialMask();
    fillFieldFromCharge();
    const std::vector<Vector_t> gathered = referenceGather();
    EXPECT_EQ(countNonFinite(gathered), 0u);
    std::vector<Vector_t> expected = gathered;
    size_t nRejected = 0;
    for (size_t i = 0; i < expected.size(); ++i) {
        if (!accept(i)) {
            expected[i] = sentinel_m;
            ++nRejected;
        }
    }
    EXPECT_GT(globalSum(nRejected), nInside / 4);

    fillFieldFromCharge();
    setEf(sentinel_m);
    bunch_m.gatherMasked(accept);
    const std::vector<Vector_t> masked = localEf();
    EXPECT_EQ(countDifferent(masked, expected), 0u);
    // and the comparison does see the rejected particles
    EXPECT_GT(countDifferent(masked, gathered), 0u);
}
