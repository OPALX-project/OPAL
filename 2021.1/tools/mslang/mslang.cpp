#include "Algorithms/Vektor.h"
#include "Utilities/MSLang.h"
#include "Utilities/MSLang/QuadTree.h"

#include "Ippl.h"
#include "Utility/IpplTimings.h"

#include <gsl/gsl_rng.h>

#include <boost/regex.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <streambuf>
#include <cstdlib>
#include <cmath>
#include <set>

Ippl *ippl;
int main(int argc, char *argv[])
{
    ippl = new Ippl(argc, argv);

    if (argc < 4) {
        std::cout << "please provide the name of the file that contains your code, the method and the number of particles" << std::endl;
        return 1;
    }
    mslang::Function *fun;

    std::ifstream t(argv[1]);
    std::string str((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());

    const unsigned int method = (atoi(argv[2]) == 0 ? 0: 1);
    const unsigned int N = atoi(argv[3]);

    // std::string str("repeat( translate(union(rectangle(0.1, 0.1), ellipse(0.1, 0.1)), -0.01, -0.02), 2, 0.1, 0.2)");

    str = boost::regex_replace(str, boost::regex("//.*?\\n"), std::string(""), boost::match_default | boost::format_all);
    str = boost::regex_replace(str, boost::regex("\\s"), std::string(""), boost::match_default | boost::format_all);

    if (parse(str, fun)) {
        // fun->print(0);
        // std::cout << "\n" << std::endl;

        std::vector<std::shared_ptr<mslang::Base> > baseBlocks;
        fun->apply(baseBlocks);

        std::ofstream out("test.gpl");
        for (std::shared_ptr<mslang::Base> & bfun: baseBlocks) {
            // bfun->print(0);
            bfun->computeBoundingBox();
            bfun->writeGnuplot(out);
        }
        out.close();

        if (baseBlocks.size() > 0) {
            Vector_t llc, urc;
            std::shared_ptr<mslang::Base> & first = baseBlocks.front();
            const mslang::BoundingBox &bb = first->bb_m;
            llc = Vector_t(bb.center_m[0] - 0.5 * bb.width_m,
                           bb.center_m[1] - 0.5 * bb.height_m,
                           0.0);
            urc = Vector_t(bb.center_m[0] + 0.5 * bb.width_m,
                           bb.center_m[1] + 0.5 * bb.height_m,
                           0.0);


            for (unsigned int i = 1; i < baseBlocks.size(); ++ i) {
                const mslang::BoundingBox &bb = baseBlocks[i]->bb_m;
                llc[0] = std::min(llc[0], bb.center_m[0] - 0.5 * bb.width_m);
                llc[1] = std::min(llc[1], bb.center_m[1] - 0.5 * bb.height_m);
                urc[0] = std::max(urc[0], bb.center_m[0] + 0.5 * bb.width_m);
                urc[1] = std::max(urc[1], bb.center_m[1] + 0.5 * bb.height_m);
            }

            double width = urc[0] - llc[0];
            double height = urc[1] - llc[1];
            llc[0] -= 1e-3 * width;
            urc[0] += 1e-3 * width;
            llc[1] -= 1e-3 * height;
            urc[1] += 1e-3 * width;

            mslang::QuadTree tree;
            tree.bb_m = mslang::BoundingBox(llc, urc);

            tree.objects_m.insert(tree.objects_m.end(), baseBlocks.begin(), baseBlocks.end());
            tree.buildUp();

            out.open("quadtree.gpl");
            tree.writeGnuplot(out);

            out.close();

            if (true) {
                out.open("particles.gpl");
                gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

                unsigned int n = 0;
                IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");
                IpplTimings::startTimer(mainTimer);

                for (unsigned int i = 0; i < N; ++ i) {
                    Vector_t X(0.0);
                    X[0] = llc[0] + (urc[0] - llc[0]) * gsl_rng_uniform(rng);
                    X[1] = llc[1] + (urc[1] - llc[1]) * gsl_rng_uniform(rng);


                    if (method == 0) {
                        if (tree.isInside(X)) {
                            ++ n;
                            out << std::setw(14) << X[0]
                                << std::setw(14) << X[1]
                                << std::endl;
                        }
                    } else {
                        for (std::shared_ptr<mslang::Base> & func: baseBlocks) {
                            if (func->isInside(X)) {
                                ++ n;
                                out << std::setw(14) << X[0]
                                    << std::setw(14) << X[1]
                                    << std::endl;
                                break;
                            }
                        }
                    }
                }
                IpplTimings::stopTimer(mainTimer);

                std::cout << (double)n / N * 100 << " % of particles passed" << std::endl;

                std::map<std::string, unsigned int> characteristicValues;
                characteristicValues.insert(std::make_pair("method", method));
                characteristicValues.insert(std::make_pair("num particles", N));
                characteristicValues.insert(std::make_pair("num base functions", baseBlocks.size()));
                std::stringstream ss;
                ss << "timing__m_=_" << method << "__np_=_" << N << "__nbf_=_" << baseBlocks.size() << ".dat";
                IpplTimings::print(ss.str(), characteristicValues);
                gsl_rng_free(rng);
            }

            for (std::shared_ptr<mslang::Base> & func: baseBlocks) {
                func.reset();
            }
        }

        out.close();
    }

    delete fun;

    return 0;
}