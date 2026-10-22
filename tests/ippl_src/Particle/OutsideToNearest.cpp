// Tests of ParticleSpatialLayout::setOutsideToNearest(), which lets update()
// give a particle outside the global domain to the node owning the point of
// the domain nearest to it, and of BinaryRepartition with such particles.
//
// The tests pass on any number of nodes, but update() swaps particles only on
// more than one node, so run them with, for example,
//   mpirun -np 4 opal_unit_tests --gtest_filter='ParticleSpatialLayout.*:ParticleBalancer.*'
#include "gtest/gtest.h"

#include "Field/Field.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Meshes/UniformCartesian.h"
#include "Particle/IntNGP.h"
#include "Particle/IpplParticleBase.h"
#include "Particle/ParticleBalancer.h"
#include "Particle/ParticleSpatialLayout.h"
#include "Utility/IpplException.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
    typedef UniformCartesian<3, double>              Mesh_t;
    typedef CenteredFieldLayout<3, Mesh_t, Cell>     FieldLayout_t;
    typedef ParticleSpatialLayout<double, 3, Mesh_t> Layout_t;
    typedef Vektor<double, 3>                        Vector_t;
    typedef NDRegion<double, 3>                      Region_t;

    const unsigned int nCells = 16;

    // Mesh and field layout, built as in FieldSolver::initCartesianFields.
    struct Grid {
        Grid():
            mesh(Index(nCells + 1), Index(nCells + 1), Index(nCells + 1)),
            layout(mesh, decomp)
        { }

        e_dim_tag decomp[3] = {PARALLEL, PARALLEL, PARALLEL};
        Mesh_t mesh;
        FieldLayout_t layout;
    };

    // Particles set up as OPAL's PartBunch: the layout is bound to the field
    // layout and makes its own mesh, whose spacing and origin are then set as
    // in PartBunch::updateFields. With the default origin, the domain is
    // [-0.85, 0.75) x [-1.7, 1.5) x [-0.4, 0.4).
    class Bunch: public IpplParticleBase<Layout_t> {
    public:
        Bunch(FieldLayout_t& fl, const Vector_t& origin = Vector_t({-0.85, -1.7, -0.4})):
            IpplParticleBase<Layout_t>(new Layout_t())
        {
            addAttribute(tag);
            addAttribute(canSwap);
            getLayout().getLayout().changeDomain(fl);
            double spacing[3] = {0.1, 0.2, 0.05};
            getMesh().set_meshSpacing(spacing);
            getMesh().set_origin(origin);
        }

        Mesh_t& getMesh() { return getLayout().getLayout().getMesh(); }
        const Region_t& getDomain() const { return getLayout().getLayout().getDomain(); }

        ParticleAttrib<int>  tag;       // index of the particle in the table of positions
        ParticleAttrib<char> canSwap;
    };

    bool isInside(const Region_t& dom, const Vector_t& x) {
        bool inside = true;
        for (unsigned int d = 0; d < 3; ++d)
            inside = inside && x[d] >= dom[d].min() && x[d] < dom[d].max();
        return inside;
    }

    bool isBitwiseEqual(const Vector_t& a, const Vector_t& b) {
        for (unsigned int d = 0; d < 3; ++d) {
            const double x = a[d], y = b[d];
            if (std::memcmp(&x, &y, sizeof(double)) != 0)
                return false;
        }
        return true;
    }

    size_t countOutside(const std::vector<Vector_t>& table, const Region_t& dom) {
        size_t n = 0;
        for (const Vector_t& x: table)
            n += !isInside(dom, x);
        return n;
    }

    // quasi-random numbers in (0, 1), the same on every node
    double quasiRandom(unsigned int k, unsigned int d) {
        const double alpha[3] = {0.6180339887498949, 0.4142135623730951, 0.7320508075688772};
        const double u = (k + 1) * alpha[d];
        return u - std::floor(u);
    }

    // The positions of the test particles. First nInside particles inside the
    // domain, denser towards its lower corner so that a repartition moves the
    // cuts. Then, if withOutside, particles outside the domain on every face
    // (from just outside, and exactly on the upper faces, which lie outside the
    // half-open domain, to 1e6 domain widths away), edge and corner.
    std::vector<Vector_t> makePositions(const Region_t& dom, bool withOutside) {
        const unsigned int nInside = 400;
        std::vector<Vector_t> table;
        Vector_t lo, width;
        for (unsigned int d = 0; d < 3; ++d) {
            lo[d] = dom[d].min();
            width[d] = dom[d].length();
        }
        for (unsigned int k = 0; k < nInside; ++k) {
            Vector_t x;
            for (unsigned int d = 0; d < 3; ++d)
                x[d] = lo[d] + width[d] * std::pow(quasiRandom(k, d), 1.5);
            table.push_back(x);
        }
        if (!withOutside)
            return table;

        unsigned int k = nInside;
        const double dist[] = {0.0, 1e-9, 0.01, 0.1, 0.5, 2.0, 1e6};  // in domain widths
        for (unsigned int d = 0; d < 3; ++d) {
            for (double s: dist) {
                Vector_t below, above;
                for (unsigned int j = 0; j < 3; ++j) {
                    below[j] = lo[j] + width[j] * quasiRandom(k++, j);
                    above[j] = lo[j] + width[j] * quasiRandom(k++, j);
                }
                below[d] = (s == 0.0 ?
                            std::nextafter(dom[d].min(), -std::numeric_limits<double>::infinity()) :
                            dom[d].min() - s * width[d]);
                above[d] = dom[d].max() + s * width[d];
                table.push_back(below);
                table.push_back(above);
            }
        }
        // edges and corners: outside in two or three dimensions
        for (unsigned int c = 0; c < 27; ++c) {
            const int side[3] = {int(c % 3) - 1, int(c / 3 % 3) - 1, int(c / 9) - 1};
            if (std::abs(side[0]) + std::abs(side[1]) + std::abs(side[2]) < 2)
                continue;
            Vector_t x;
            for (unsigned int j = 0; j < 3; ++j) {
                x[j] = lo[j] + width[j] * quasiRandom(k++, j);
                if (side[j] < 0)
                    x[j] = dom[j].min() - 0.3 * width[j];
                else if (side[j] > 0)
                    x[j] = dom[j].max() + 0.7 * width[j];
            }
            table.push_back(x);
        }
        return table;
    }

    // Append n particles outside the domain, each beyond one of the upper
    // faces. Unlike those of makePositions(), they lie on one side of the
    // domain only, so counting them in a balance at the nearest points of the
    // domain would move the cuts towards the upper faces.
    void addOneSided(std::vector<Vector_t>& table, const Region_t& dom, unsigned int n) {
        const unsigned int k0 = table.size();
        for (unsigned int k = k0; k < k0 + n; ++k) {
            Vector_t x;
            for (unsigned int j = 0; j < 3; ++j)
                x[j] = dom[j].min() + dom[j].length() * quasiRandom(k, j);
            const unsigned int d = k % 3;
            x[d] = dom[d].max() + 0.3 * dom[d].length();
            table.push_back(x);
        }
    }

    // the table mirrored through the centre of the domain
    std::vector<Vector_t> mirror(const std::vector<Vector_t>& table, const Region_t& dom) {
        std::vector<Vector_t> mirrored(table);
        for (Vector_t& x: mirrored)
            for (unsigned int d = 0; d < 3; ++d)
                x[d] = (dom[d].min() + dom[d].max()) - x[d];
        return mirrored;
    }

    // create the particles of the table on node 0, tagged with their index
    void fill(Bunch& bunch, const std::vector<Vector_t>& table) {
        if (Ippl::myNode() != 0)
            return;
        bunch.create(table.size());
        for (size_t k = 0; k < table.size(); ++k) {
            bunch.R[k] = table[k];
            bunch.tag[k] = k;
            bunch.canSwap[k] = 1;
        }
    }

    // Is x on this node as update() should place it? A particle inside the
    // domain must lie in a local rnode, tested half-open as update() does. A
    // particle outside it must be on a node whose rnode touches the point of
    // the domain nearest to it; that point lies on the domain boundary, so the
    // closed rnode is tested, with a tolerance far below one cell.
    bool ownsNearestPoint(Layout_t& layout, const Region_t& dom, const Vector_t& x) {
        const bool inside = isInside(dom, x);
        Vector_t q;
        for (unsigned int d = 0; d < 3; ++d)
            q[d] = std::min(std::max(x[d], dom[d].min()), dom[d].max());

        RegionLayout<double, 3, Mesh_t>& rl = layout.getLayout();
        for (auto it = rl.begin_iv(); it != rl.end_iv(); ++it) {
            const Region_t& r = (*it).second->getDomain();
            bool owns = true;
            for (unsigned int d = 0; d < 3; ++d) {
                const double tol = 1e-9 * dom[d].length();
                if (inside)
                    owns = owns && q[d] >= r[d].min() && q[d] < r[d].max();
                else
                    owns = owns && q[d] >= r[d].min() - tol && q[d] <= r[d].max() + tol;
            }
            if (owns)
                return true;
        }
        return false;
    }

    // Number of failures, summed over all nodes, of the checks after an
    // update: every particle of the table exists exactly once on all nodes,
    // its position is bitwise the one of the table, and it is on the node
    // that ownsNearestPoint() expects. Particles whose entry in pinned is set
    // must instead still be on node 0.
    int checkLayout(Bunch& bunch, const std::vector<Vector_t>& table,
                    const std::vector<char>& pinned = std::vector<char>()) {
        const Layout_t::ParticlePos_t& R = bunch.R;
        const ParticleAttrib<int>& tag = bunch.tag;
        std::vector<int> seen(table.size(), 0);
        int failures = 0;
        for (size_t i = 0; i < bunch.getLocalNum(); ++i) {
            const int k = tag[i];
            if (k < 0 || k >= (int)table.size()) {
                ++failures;
                continue;
            }
            ++seen[k];
            if (!isBitwiseEqual(R[i], table[k]))
                ++failures;
            if (!pinned.empty() && pinned[k]) {
                failures += (Ippl::myNode() != 0);
            } else if (!ownsNearestPoint(bunch.getLayout(), bunch.getDomain(), R[i])) {
                ++failures;
            }
        }
        allreduce(&failures, 1, std::plus<int>());
        allreduce(seen.data(), seen.size(), std::plus<int>());
        for (int s: seen)
            failures += (s != 1);
        return failures;
    }

    // the outside counts of all nodes, summed
    size_t getGlobalOutsideCount(Bunch& bunch) {
        size_t count = bunch.getLayout().getOutsideCount();
        allreduce(&count, 1, std::plus<size_t>());
        return count;
    }

    // the outside count update() should report; without a swap it is 0
    size_t expectedOutsideCount(size_t nOutside) {
        return (Ippl::getNodes() > 1 ? nOutside : 0);
    }

    std::vector<int> getLocalTags(Bunch& bunch) {
        std::vector<int> tags(bunch.tag.begin(), bunch.tag.end());
        std::sort(tags.begin(), tags.end());
        return tags;
    }

    // the local vnodes of a field layout, as sorted strings
    std::vector<std::string> getLocalVnodes(FieldLayout_t& fl) {
        std::vector<std::string> vnodes;
        for (auto it = fl.begin_iv(); it != fl.end_iv(); ++it) {
            std::ostringstream domain;
            domain << (*it).second->getDomain();
            vnodes.push_back(domain.str());
        }
        std::sort(vnodes.begin(), vnodes.end());
        return vnodes;
    }
}

TEST(ParticleSpatialLayout, OutsideThrowsByDefault)
{
    Grid grid;
    Bunch bunch(grid.layout);
    EXPECT_FALSE(bunch.getLayout().getOutsideToNearest());

    const std::vector<Vector_t> table = makePositions(bunch.getDomain(), true);
    fill(bunch, table);
    if (Ippl::getNodes() > 1) {
        EXPECT_THROW(bunch.update(), IpplException);
    } else {
        // one node does not swap and does not look at the positions at all
        EXPECT_NO_THROW(bunch.update());
    }
}

TEST(ParticleSpatialLayout, OutsideToNearest)
{
    Grid grid;
    Bunch bunch(grid.layout);
    const Region_t& dom = bunch.getDomain();
    const std::vector<Vector_t> table = makePositions(dom, true);
    const size_t nOutside = countOutside(table, dom);
    EXPECT_EQ(nOutside, 3u * 14u + 20u);

    bunch.getLayout().setOutsideToNearest(true);
    EXPECT_TRUE(bunch.getLayout().getOutsideToNearest());

    // from node 0 to the nodes owning the particles or the nearest points
    fill(bunch, table);
    ASSERT_NO_THROW(bunch.update());
    EXPECT_EQ(bunch.getTotalNum(), table.size());
    EXPECT_EQ(checkLayout(bunch, table), 0);
    EXPECT_EQ(getGlobalOutsideCount(bunch), expectedOutsideCount(nOutside));

    // nothing moves: every particle stays where it is
    const std::vector<int> tags = getLocalTags(bunch);
    ASSERT_NO_THROW(bunch.update());
    EXPECT_EQ(getLocalTags(bunch), tags);
    EXPECT_EQ(checkLayout(bunch, table), 0);
    EXPECT_EQ(getGlobalOutsideCount(bunch), expectedOutsideCount(nOutside));

    // the outside particles move to the opposite faces, and to other nodes
    const std::vector<Vector_t> mirrored = mirror(table, dom);
    for (size_t i = 0; i < bunch.getLocalNum(); ++i)
        bunch.R[i] = mirrored[bunch.tag[i]];
    ASSERT_NO_THROW(bunch.update());
    EXPECT_EQ(bunch.getTotalNum(), table.size());
    EXPECT_EQ(checkLayout(bunch, mirrored), 0);
    EXPECT_EQ(getGlobalOutsideCount(bunch),
              expectedOutsideCount(countOutside(mirrored, dom)));

    // and without the policy the same particles are refused again
    bunch.getLayout().setOutsideToNearest(false);
    if (Ippl::getNodes() > 1) {
        EXPECT_THROW(bunch.update(), IpplException);
    }
}

TEST(ParticleSpatialLayout, OutsideToNearestCanSwap)
{
    Grid grid;
    Bunch bunch(grid.layout);
    const Region_t& dom = bunch.getDomain();
    const std::vector<Vector_t> table = makePositions(dom, true);

    // outside particles that may not be swapped are not looked at
    fill(bunch, table);
    std::vector<char> pinned(table.size(), 0);
    for (size_t k = 0; k < table.size(); ++k)
        pinned[k] = !isInside(dom, table[k]);
    for (size_t i = 0; i < bunch.getLocalNum(); ++i)
        bunch.canSwap[i] = !pinned[bunch.tag[i]];
    ASSERT_NO_THROW(bunch.update(bunch.canSwap));
    EXPECT_EQ(checkLayout(bunch, table, pinned), 0);

    // let every other one of them be swapped: that is refused without the
    // policy and gives them to the nodes owning the nearest points with it
    size_t nSwapped = 0;
    for (size_t k = 0; k < table.size(); ++k) {
        if (pinned[k] && k % 2 == 0) {
            pinned[k] = 0;
            ++nSwapped;
        }
    }
    for (size_t i = 0; i < bunch.getLocalNum(); ++i)
        bunch.canSwap[i] = !pinned[bunch.tag[i]];
    if (Ippl::getNodes() > 1) {
        EXPECT_THROW(bunch.update(bunch.canSwap), IpplException);
    }

    bunch.getLayout().setOutsideToNearest(true);
    ASSERT_NO_THROW(bunch.update(bunch.canSwap));
    EXPECT_EQ(bunch.getTotalNum(), table.size());
    EXPECT_EQ(checkLayout(bunch, table, pinned), 0);
    EXPECT_EQ(getGlobalOutsideCount(bunch), expectedOutsideCount(nSwapped));
}

TEST(ParticleSpatialLayout, OutsideToNearestNotFinite)
{
    // a particle with a coordinate that is NaN or infinite is still refused,
    // also when another of its coordinates lies outside the domain
    const double values[] = {std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity()};
    for (double v: values) {
        for (unsigned int d = 0; d < 3; ++d) {
            Grid grid;
            Bunch bunch(grid.layout);
            const Region_t& dom = bunch.getDomain();
            bunch.getLayout().setOutsideToNearest(true);
            std::vector<Vector_t> table = makePositions(dom, false);
            const unsigned int e = (d + 1) % 3;
            table.back()[d] = v;
            table.back()[e] = dom[e].max() + 0.5 * dom[e].length();
            fill(bunch, table);
            if (Ippl::getNodes() > 1) {
                EXPECT_THROW(bunch.update(), IpplException) << "coordinate " << d << " = " << v;
            } else {
                EXPECT_NO_THROW(bunch.update());
            }
        }
    }
}

TEST(ParticleSpatialLayout, OutsideToNearestFarFromOrigin)
{
    // A domain about 1e10 of its widths away from the origin. There, 1e-7 of
    // the width is below the resolution of the coordinates, so a shift by that
    // alone would leave the lookup point of a particle beyond an upper face on
    // the face, outside the half-open domain.
    Grid grid;
    Bunch bunch(grid.layout, Vector_t({1e10, -2e10, 5e10}));
    const Region_t& dom = bunch.getDomain();
    for (unsigned int d = 0; d < 3; ++d)
        ASSERT_EQ(dom[d].max() - 1e-7 * dom[d].length(), dom[d].max());

    const std::vector<Vector_t> table = makePositions(dom, true);
    fill(bunch, table);
    if (Ippl::getNodes() > 1) {
        EXPECT_THROW(bunch.update(), IpplException);
    }

    bunch.getLayout().setOutsideToNearest(true);
    ASSERT_NO_THROW(bunch.update());
    EXPECT_EQ(bunch.getTotalNum(), table.size());
    EXPECT_EQ(checkLayout(bunch, table), 0);
    EXPECT_EQ(getGlobalOutsideCount(bunch), expectedOutsideCount(countOutside(table, dom)));
}

TEST(ParticleSpatialLayout, OutsideToNearestNoCaching)
{
    // the ghost particles are selected by the rnodes, so caching is refused
    // together with the policy, in either order
    Grid grid;
    Bunch bunch(grid.layout);
    Layout_t& layout = bunch.getLayout();
    layout.enableCaching();
    EXPECT_THROW(layout.setOutsideToNearest(true), IpplException);
    EXPECT_FALSE(layout.getOutsideToNearest());
    EXPECT_NO_THROW(layout.setOutsideToNearest(false));
    layout.disableCaching();

    layout.setOutsideToNearest(true);
    EXPECT_THROW(layout.enableCaching(), IpplException);
    layout.setOutsideToNearest(false);
    EXPECT_NO_THROW(layout.enableCaching());
    layout.disableCaching();
}

TEST(ParticleBalancer, BinaryRepartitionOutside)
{
    // two identical grids: one with all particles and the policy, one with
    // only the particles inside the domain and without it. The outside
    // particles are those of makePositions() and, beyond the upper faces, as
    // many more as there are inside.
    Grid gridAll, gridInside;
    Bunch all(gridAll.layout), inside(gridInside.layout);
    const Region_t& dom = all.getDomain();
    const std::vector<Vector_t> tableInside = makePositions(dom, false);
    std::vector<Vector_t> table = makePositions(dom, true);
    addOneSided(table, dom, tableInside.size());
    const size_t nOutside = countOutside(table, dom);
    ASSERT_EQ(table.size(), tableInside.size() + nOutside);

    all.getLayout().setOutsideToNearest(true);
    fill(all, table);
    fill(inside, tableInside);
    ASSERT_NO_THROW(all.update());
    ASSERT_NO_THROW(inside.update());
    const std::vector<std::string> vnodesBefore = getLocalVnodes(gridAll.layout);

    EXPECT_TRUE(BinaryRepartition(all));
    EXPECT_TRUE(BinaryRepartition(inside));

    // the particles outside the domain do not count in the balance
    const std::vector<std::string> vnodes = getLocalVnodes(gridAll.layout);
    EXPECT_EQ(vnodes, getLocalVnodes(gridInside.layout));
    int changed = (vnodes != vnodesBefore);
    allreduce(&changed, 1, std::plus<int>());
    EXPECT_EQ(changed > 0, Ippl::getNodes() > 1);

    // the stencil of IPPL's gradient needs vnodes of at least two cells
    for (auto it = gridAll.layout.begin_iv(); it != gridAll.layout.end_iv(); ++it)
        for (unsigned int d = 0; d < 3; ++d)
            EXPECT_GE((*it).second->getDomain()[d].length(), 2u);

    // no particle is lost or moved, and each is on its new node
    EXPECT_EQ(all.getTotalNum(), table.size());
    EXPECT_EQ(checkLayout(all, table), 0);
    EXPECT_EQ(getGlobalOutsideCount(all), expectedOutsideCount(nOutside));
    EXPECT_EQ(checkLayout(inside, tableInside), 0);

    // the number density that BinaryRepartition balances counts exactly the
    // particles inside the domain
    Field<double, 3, Mesh_t, Cell> density(all.getMesh(), gridAll.layout, GuardCellSizes<3>(1));
    Field<double, 3, Mesh_t, Cell> densityInside(inside.getMesh(), gridInside.layout,
                                                 GuardCellSizes<3>(1));
    density = 0.0;
    densityInside = 0.0;
    scatterInsideDomain(density, all.R, dom, IntNGP());
    scatter(densityInside, inside.R, IntNGP());
    EXPECT_EQ(sum(density), double(tableInside.size()));
    EXPECT_EQ(sum(density * density), sum(densityInside * densityInside));
}

TEST(ParticleBalancer, BinaryRepartitionOutsideThinVnode)
{
    // The balancer refuses a partition with a vnode one cell wide, whatever
    // the weights. With every inside particle in the lowest cell, the first
    // cut would leave one, so the repartition fails and changes nothing; the
    // outside particles do not change that.
    Grid grid;
    Bunch bunch(grid.layout);
    const Region_t& dom = bunch.getDomain();
    std::vector<Vector_t> table = makePositions(dom, true);
    for (Vector_t& x: table) {
        if (!isInside(dom, x))
            continue;
        for (unsigned int d = 0; d < 3; ++d)
            x[d] = dom[d].min() + (x[d] - dom[d].min()) / nCells;
    }
    const size_t nOutside = countOutside(table, dom);

    bunch.getLayout().setOutsideToNearest(true);
    fill(bunch, table);
    ASSERT_NO_THROW(bunch.update());
    const std::vector<std::string> vnodesBefore = getLocalVnodes(grid.layout);
    const std::vector<int> tags = getLocalTags(bunch);

    EXPECT_EQ(BinaryRepartition(bunch), Ippl::getNodes() == 1);
    EXPECT_EQ(getLocalVnodes(grid.layout), vnodesBefore);
    EXPECT_EQ(getLocalTags(bunch), tags);
    EXPECT_EQ(bunch.getTotalNum(), table.size());
    EXPECT_EQ(checkLayout(bunch, table), 0);
    EXPECT_EQ(getGlobalOutsideCount(bunch), expectedOutsideCount(nOutside));
}
