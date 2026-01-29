#ifndef MSLANG_MASK_H
#define MSLANG_MASK_H

#include "Utilities/MSLang.h"
#include "Utilities/MSLang/Rectangle.h"

namespace mslang {
    struct Mask: public Function {
        static bool parse_detail(iterator &it, const iterator &end, Function* &fun);
        virtual void print(int ident);
        virtual void apply(std::vector<std::shared_ptr<Base> > &bfuncs);

        std::vector<Rectangle> pixels_m;

    private:

        struct IntPoint {
            IntPoint(int x, int y):
                x_m(x),
                y_m(y)
            { }

            int x_m;
            int y_m;
        };

        typedef std::pair<IntPoint, IntPoint> IntPixel_t;

        std::vector<IntPixel_t> minimizeNumberOfRectangles(std::vector<bool> pixels,
                                                           unsigned int height,
                                                           unsigned int width);

        std::pair<IntPoint, IntPoint> findMaximalRectangle(const std::vector<bool> &pixels,
                                                           unsigned int height,
                                                           unsigned int width) const;

        unsigned int computeArea(const IntPoint &ll, const IntPoint &ur) const;

        void updateCache(const std::vector<bool> &pixels, std::vector<unsigned int> &cache, unsigned int y) const;
    };
}

#endif