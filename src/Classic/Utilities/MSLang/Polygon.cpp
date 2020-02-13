#include "Utilities/MSLang/Polygon.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"
#include "Utilities/Mesher.h"
#include "Physics/Physics.h"

#include <boost/regex.hpp>

namespace mslang {
    void Polygon::triangulize(std::vector<Vector_t> &nodes) {
        Mesher mesher(nodes);
        triangles_m = mesher.getTriangles();
    }

    bool Polygon::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Polygon *poly = static_cast<Polygon*>(fun);

        ArgumentExtractor arguments(std::string(it, end));
        std::vector<Vector_t> nodes;

        for (unsigned int i = 0; i + 1 < arguments.getNumArguments(); i += 2) {
            try {
                double x = parseMathExpression(arguments.get(i));
                double y = parseMathExpression(arguments.get(i + 1));
                nodes.push_back(Vector_t(x, y, 1.0));
            } catch (std::runtime_error &e) {
                std::cout << e.what() << std::endl;
                return false;
            }
        }

        if (nodes.size() < 3) return false;

        poly->triangulize(nodes);

        it += (arguments.getLengthConsumed() + 1);
        return true;
    }

    void Polygon::print(int /*ident*/) {
        // for (auto pix: pixels_m) pix.print(ident);
    }

    void Polygon::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        for (Triangle &tri: triangles_m)
            bfuncs.push_back(std::shared_ptr<Base>(tri.clone()));
    }
}
