#include "Utilities/MSLang.h"
#include "Utilities/MSLang/Ellipse.h"
#include "Utilities/MSLang/Rectangle.h"
#include "Utilities/MSLang/Triangle.h"
#include "Utilities/MSLang/Polygon.h"
#include "Utilities/MSLang/Mask.h"
#include "Utilities/MSLang/Translation.h"
#include "Utilities/MSLang/Rotation.h"
#include "Utilities/MSLang/Shear.h"
#include "Utilities/MSLang/Repeat.h"
#include "Utilities/MSLang/Union.h"
#include "Utilities/MSLang/Difference.h"
#include "Utilities/MSLang/SymmetricDifference.h"
#include "Utilities/MSLang/Intersection.h"
#include "Utilities/MSLang/QuadTree.h"
#include "Utilities/PortableBitmapReader.h"
#include "Utilities/Mesher.h"
#include "Algorithms/Quaternion.h"
#include "Physics/Physics.h"


#include <boost/regex.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <streambuf>
#include <cstdlib>
#include <cmath>

namespace mslang {
    const std::string Function::UDouble  = "([0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)";
    const std::string Function::Double = "(-?[0-9]+\\.?[0-9]*([Ee][+-]?[0-9]+)?)";
    const std::string Function::UInt = "([0-9]+)";
    const std::string Function::FCall = "([a-z_]*)\\((.*)";

    bool parse(std::string str, Function* &fun) {
        iterator it = str.begin();
        iterator end = str.end();
        if (!Function::parse(it, end, fun)) {
            std::cout << "parsing failed here:" << std::string(it, end) << std::endl;
            return false;
        }

        return true;
    }

    bool Function::parse(iterator &it, const iterator &end, Function* &fun) {
        boost::regex functionCall(Function::FCall);
        boost::smatch what;

        std::string str(it, end);
        if( !boost::regex_match(str , what, functionCall ) ) return false;

        std::string identifier = what[1];
        std::string arguments = what[2];
        unsigned int shift = identifier.size() + 1;

        if (identifier == "rectangle") {
            fun = new Rectangle;
            it += shift;
            if (!Rectangle::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "ellipse") {
            fun = new Ellipse;
            it += shift;
            if (!Ellipse::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "polygon") {
            fun = new Polygon;
            it += shift;
            if (!Polygon::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "mask") {
            fun = new Mask;
            it += shift;

            return Mask::parse_detail(it, end, fun);
        } else if (identifier == "repeat") {
            fun = new Repeat;
            it += shift;
            if (!Repeat::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "rotate") {
            fun = new Rotation;
            it += shift;
            if (!Rotation::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "translate") {
            fun = new Translation;
            it += shift;
            if (!Translation::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "shear") {
            fun = new Shear;
            it += shift;
            if (!Shear::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "union") {
            fun = new Union;
            it += shift;
            if (!Union::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "difference") {
            fun = new Difference;
            it += shift;
            if (!Difference::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "symmetric_difference") {

            fun = new SymmetricDifference;
            it += shift;
            if (!SymmetricDifference::parse_detail(it, end, fun)) return false;

            return true;
        } else if (identifier == "intersection") {
            fun = new Intersection;
            it += shift;
            if (!Intersection::parse_detail(it, end, fun)) return false;

            return true;
        }

        return (it == end);
    }
}