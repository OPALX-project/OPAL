//
// Tests for FFTPoissonSolver (open boundaries, INTEGRATED and STANDARD Green's function)
//
// The solver returns phi_i = sum_j rho_j G_ij, with no factor of the cell volume
// and no 1/(4 pi eps0). A unit charge in one cell therefore returns the Green's
// function itself, which is checked against independent references:
//  - INTEGRATED: the average of 1/r over a source cell (Qiang et al.), evaluated
//    here by midpoint sub-sampling. It is finite in the self cell (2.38008/h for a
//    cube), so no neighbour value may be substituted for it.
//  - STANDARD: the point kernel 1/r, regularised at r = 0 by the neighbour value.
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

#include "Solvers/FFTPoissonSolver.h"

#include "opal_test_utilities/SilenceTest.h"

#include <array>
#include <cmath>
#include <functional>
#include <memory>
#include <string>

namespace {
    constexpr int nCells = 12;
    constexpr int centre = nCells / 2;

    using Offset = std::array<int, 3>;

    // Mean of 1/|v + u| over u in the cell [-h/2, h/2]^3, v = (i hx, j hy, k hz).
    double cellAverage(const Offset& offset, const Vector_t& h, int nSub = 60) {
        double sum = 0.0;
        for (int a = 0; a < nSub; ++a) {
            const double x = offset[0] * h(0) + ((a + 0.5) / nSub - 0.5) * h(0);
            for (int b = 0; b < nSub; ++b) {
                const double y = offset[1] * h(1) + ((b + 0.5) / nSub - 0.5) * h(1);
                for (int c = 0; c < nSub; ++c) {
                    const double z = offset[2] * h(2) + ((c + 0.5) / nSub - 0.5) * h(2);
                    sum += 1.0 / std::sqrt(x * x + y * y + z * z);
                }
            }
        }
        return sum / (double(nSub) * nSub * nSub);
    }

    double pointKernel(const Offset& offset, const Vector_t& h) {
        const double x = offset[0] * h(0);
        const double y = offset[1] * h(1);
        const double z = offset[2] * h(2);
        return 1.0 / std::sqrt(x * x + y * y + z * z);
    }

    class FFTPoissonSolverTest: public ::testing::Test {
    public:
        FFTPoissonSolverTest():
            mesh_m(meshDomain()),
            layout_m(mesh_m, decomp_m)
        { }

        void SetUp() override {
            setSpacing(Vector_t({1.0, 1.0, 1.0}));
            for (int i = 0; i < 2 * 3; ++i) {
                bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
            }
            rho_m.initialize(mesh_m, layout_m, GuardCellSizes<3>(1), bc_m);
        }

    protected:
        static NDIndex<3> meshDomain() {
            NDIndex<3> domain;
            for (unsigned d = 0; d < 3; ++d) {
                domain[d] = Index(nCells + 1);
            }
            return domain;
        }

        void setSpacing(const Vector_t& h) {
            hr_m = h;
            mesh_m.set_meshSpacing(&(hr_m[0]));
            mesh_m.set_origin(Vector_t({0.0, 0.0, 0.0}));
        }

        std::unique_ptr<FFTPoissonSolver> makeSolver(const std::string& greens) {
            return std::unique_ptr<FFTPoissonSolver>(
                new FFTPoissonSolver(&mesh_m, &layout_m, greens, "OPEN"));
        }

        static NDIndex<3> point(const Offset& c) {
            NDIndex<3> p;
            for (unsigned d = 0; d < 3; ++d) {
                p[d] = Index(c[d], c[d]);
            }
            return p;
        }

        void clear() {
            rho_m = 0.0;
        }

        void setCharge(const Offset& c, double q) {
            const NDIndex<3> p = point(c);
            if (layout_m.getLocalNDIndex().contains(p)) {
                rho_m.localElement(p) = q;
            }
        }

        // Value of the cell on whichever rank owns it, on every rank.
        double valueAt(const Offset& c) {
            const NDIndex<3> p = point(c);
            double v = 0.0;
            if (layout_m.getLocalNDIndex().contains(p)) {
                v = rho_m.localElement(p);
            }
            allreduce(v, 1, std::plus<double>());
            return v;
        }

        // Potential at centre + offset of a unit charge at the centre cell.
        double potential(const Offset& offset) {
            return valueAt({centre + offset[0], centre + offset[1], centre + offset[2]});
        }

        void solveUnitCharge(FFTPoissonSolver& solver) {
            clear();
            setCharge({centre, centre, centre}, 1.0);
            solver.computePotential(rho_m, hr_m);
        }

        e_dim_tag decomp_m[3] = {SERIAL, SERIAL, PARALLEL};
        Mesh_t mesh_m;
        FieldLayout_t layout_m;
        BConds<double, 3, Mesh_t, Center_t> bc_m;
        Field_t rho_m;
        Vector_t hr_m;
    };

    const std::array<Offset, 8> probeCells = {{
        {{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}},
        {{1, 1, 0}}, {{1, 1, 1}}, {{2, 0, 1}}, {{-3, 2, -1}}
    }};
}

TEST_F(FFTPoissonSolverTest, IntegratedSelfCellIsCellAverage) {
    OpalTestUtilities::SilenceTest silencer;

    // A cube of side h: <1/r> = 2.38008 / h. The pre-fix code replaced it by the value
    // of the neighbour along z (0.98759 / h).
    auto solver = makeSolver("INTEGRATED");
    solveUnitCharge(*solver);

    const double self = potential({0, 0, 0});
    EXPECT_NEAR(self, 2.38008, 1e-3);
    EXPECT_GT(self, 2.0 * potential({0, 0, 1}));
}

TEST_F(FFTPoissonSolverTest, IntegratedMatchesCellAverageOfPointKernel) {
    OpalTestUtilities::SilenceTest silencer;

    auto solver = makeSolver("INTEGRATED");
    solveUnitCharge(*solver);

    for (const Offset& c : probeCells) {
        const double expected = cellAverage(c, hr_m);
        EXPECT_NEAR(potential(c), expected, 1e-3 * expected)
            << "offset " << c[0] << " " << c[1] << " " << c[2];
    }
}

TEST_F(FFTPoissonSolverTest, IntegratedAnisotropicCells) {
    OpalTestUtilities::SilenceTest silencer;

    setSpacing(Vector_t({0.5, 1.0, 2.0}));
    auto solver = makeSolver("INTEGRATED");
    solveUnitCharge(*solver);

    for (const Offset& c : probeCells) {
        const double expected = cellAverage(c, hr_m);
        EXPECT_NEAR(potential(c), expected, 1e-3 * expected)
            << "offset " << c[0] << " " << c[1] << " " << c[2];
    }
}

TEST_F(FFTPoissonSolverTest, IntegratedScalesInverselyWithSpacing) {
    OpalTestUtilities::SilenceTest silencer;

    // The solver is reused with a new spacing, as the tracker does every step.
    auto solver = makeSolver("INTEGRATED");
    solveUnitCharge(*solver);
    std::array<double, 8> reference;
    for (size_t n = 0; n < probeCells.size(); ++n) {
        reference[n] = potential(probeCells[n]);
    }

    setSpacing(Vector_t({2.0, 2.0, 2.0}));
    solveUnitCharge(*solver);
    for (size_t n = 0; n < probeCells.size(); ++n) {
        EXPECT_NEAR(potential(probeCells[n]), 0.5 * reference[n], 1e-12 * reference[n]);
    }
}

TEST_F(FFTPoissonSolverTest, IntegratedFarFieldApproachesPointCharge) {
    OpalTestUtilities::SilenceTest silencer;

    auto solver = makeSolver("INTEGRATED");
    solveUnitCharge(*solver);

    const Offset far = {5, 4, 3};
    EXPECT_NEAR(potential(far), pointKernel(far, hr_m), 5e-3 * pointKernel(far, hr_m));
}

TEST_F(FFTPoissonSolverTest, IntegratedIsSymmetricAndLinear) {
    OpalTestUtilities::SilenceTest silencer;

    auto solver = makeSolver("INTEGRATED");
    solveUnitCharge(*solver);

    // Mirror symmetry of the kernel in every axis, including the self-cell neighbours.
    for (const Offset& c : probeCells) {
        const Offset m = {-c[0], -c[1], -c[2]};
        if (centre + std::abs(c[0]) < nCells && centre + std::abs(c[1]) < nCells
            && centre + std::abs(c[2]) < nCells) {
            const double v = potential(c);
            EXPECT_NEAR(potential(m), v, 1e-12 * v);
        }
    }

    // Two charges: the potential is the sum of the two single-charge potentials.
    const Offset qa = {centre, centre, centre};
    const Offset qb = {centre + 2, centre - 1, centre + 1};
    const Offset probe = {centre - 2, centre + 1, centre - 3};
    const double single = potential({-2, 1, -3});
    clear();
    setCharge(qb, 1.0);
    solver->computePotential(rho_m, hr_m);
    const double other = valueAt(probe);
    clear();
    setCharge(qa, 1.0);
    setCharge(qb, 1.0);
    solver->computePotential(rho_m, hr_m);
    EXPECT_NEAR(valueAt(probe), single + other, 1e-12 * (single + other));
}

TEST_F(FFTPoissonSolverTest, StandardIsPointKernelWithRegularisedSelfCell) {
    OpalTestUtilities::SilenceTest silencer;

    setSpacing(Vector_t({0.5, 1.0, 2.0}));
    auto solver = makeSolver("STANDARD");
    solveUnitCharge(*solver);

    for (const Offset& c : probeCells) {
        if (c[0] == 0 && c[1] == 0 && c[2] == 0) {
            continue;
        }
        const double expected = pointKernel(c, hr_m);
        EXPECT_NEAR(potential(c), expected, 1e-9 * expected)
            << "offset " << c[0] << " " << c[1] << " " << c[2];
    }
    // The point kernel diverges at r = 0 and is replaced by the value of the z-neighbour.
    EXPECT_NEAR(potential({0, 0, 0}), potential({0, 0, 1}), 1e-12);
}
