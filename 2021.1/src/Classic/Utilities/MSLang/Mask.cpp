#include "Utilities/MSLang/Mask.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"
#include "Utilities/PortableBitmapReader.h"

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

namespace mslang {
    void Mask::updateCache(const std::vector<bool> &pixels, std::vector<unsigned int> &cache, unsigned int y) const {
        const unsigned int M = cache.size();
        unsigned int idx = y * M;
        for (unsigned int x = 0; x < M; ++ x, ++ idx) {
            if (pixels[idx]) {
                ++ cache[x];
            } else {
                cache[x] = 0;
            }
        }
    }

    unsigned int Mask::computeArea(const Mask::IntPoint &ll, const Mask::IntPoint &ur) const {
        if ((ur.x_m > ll.x_m) && (ur.y_m > ll.y_m))
            return (ur.x_m - ll.x_m) * (ur.y_m - ll.y_m);

        return 0;
    }

    std::pair<Mask::IntPoint, Mask::IntPoint> Mask::findMaximalRectangle(const std::vector<bool> &pixels,
                                                                         unsigned int N, /* height */
                                                                         unsigned int M /* width */) const {

        // This algorithm was presented in
        // http://www.drdobbs.com/database/the-maximal-rectangle-problem/184410529
        // by David Vandevoorde, April 01, 1998

        unsigned int bestArea = 0;
        IntPoint bestLL(0, 0), bestUR(0,0);
        std::vector<unsigned int> cache(M, 0);
        std::stack<std::pair<unsigned int, unsigned int> > stack;
        for (unsigned int y = N - 1; y + 1 > 0; -- y) {
            updateCache(pixels, cache, y);
            unsigned int height = 0;
            for (unsigned int x = 0; x < M; ++ x) {
                if (cache[x] > height) {
                    stack.push(std::make_pair(x, height));
                    height = cache[x];
                } else if (cache[x] < height) {
                    std::pair<unsigned int, unsigned int> tmp;
                    do {
                        tmp = stack.top();
                        stack.pop();
                        if (x > tmp.first && height * (x - tmp.first) > bestArea) {
                            bestLL.x_m = tmp.first; bestLL.y_m = y;
                            bestUR.x_m = x; bestUR.y_m = y + height;
                            bestArea = height * (x - tmp.first);
                        }
                        height = tmp.second;
                    } while (stack.size() > 0 && cache[x] < height);
                    height = cache[x];
                    if (height != 0) {
                        stack.push(std::make_pair(tmp.first, height));
                    }
                }
            }

            if (stack.size() > 0) {
                std::pair<unsigned int, unsigned int> tmp = stack.top();
                stack.pop();
                if (M > tmp.first && height * (M - tmp.first) > bestArea) {
                    bestLL.x_m = tmp.first; bestLL.y_m = y;
                    bestUR.x_m = M; bestUR.y_m = y + height;
                    bestArea = height * (M - tmp.first);
                }
            }
        }

        return std::make_pair(bestLL, bestUR);
    }

    std::vector<Mask::IntPixel_t> Mask::minimizeNumberOfRectangles(std::vector<bool> pixels,
                                                                   unsigned int N,/* height */
                                                                   unsigned int M /* width */)
    {
        std::vector<IntPixel_t> rectangles;

        unsigned int maxArea = 0;
        while (true) {
            IntPixel_t pix = findMaximalRectangle(pixels, N, M);
            unsigned int area = computeArea(pix.first, pix.second);
            if (area > maxArea) maxArea = area;
            if (1000 * area < maxArea || area <= 1) {
                break;
            }

            rectangles.push_back(pix);

            for (int y = pix.first.y_m; y < pix.second.y_m; ++ y) {
                unsigned int idx = y * M + pix.first.x_m;
                for (int x = pix.first.x_m; x < pix.second.x_m; ++ x, ++ idx) {
                    pixels[idx] = false;
                }
            }
        }

        unsigned int idx = 0;
        for (unsigned int y = 0; y < N; ++ y) {
            for (unsigned int x = 0; x < M; ++ x, ++idx) {
                if (pixels[idx]) {
                    IntPoint ll(x, y);
                    IntPoint ur(x + 1, y + 1);
                    rectangles.push_back(IntPixel_t(ll, ur));
                }
            }
        }

        return rectangles;
    }

    bool Mask::parse_detail(iterator &it, const iterator &end, Function* &fun) {
        Mask *pixmap = static_cast<Mask*>(fun);

        ArgumentExtractor arguments(std::string(it, end));
        std::string filename = arguments.get(0);
        if (filename[0] == '\'' && filename.back() == '\'') {
            filename = filename.substr(1, filename.length() - 2);
        }

        if (!boost::filesystem::exists(filename)) {
            ERRORMSG("file '" << filename << "' doesn't exists" << endl);
            return false;
        }

        PortableBitmapReader reader(filename);
        unsigned int width = reader.getWidth();
        unsigned int height = reader.getHeight();

        double pixel_width;
        double pixel_height;
        try {
            pixel_width = parseMathExpression(arguments.get(1)) / width;
            pixel_height = parseMathExpression(arguments.get(2)) / height;
        } catch (std::runtime_error &e) {
            std::cout << e.what() << std::endl;
            return false;
        }

        if (pixel_width < 0.0) {
            std::cout << "Mask: a negative width provided '"
                      << arguments.get(0) << " = " << pixel_width * width << "'"
                      << std::endl;
            return false;
        }

        if (pixel_height < 0.0) {
            std::cout << "Mask: a negative height provided '"
                      << arguments.get(1) << " = " << pixel_height * height << "'"
                      << std::endl;
            return false;
        }

        auto maxRect = pixmap->minimizeNumberOfRectangles(reader.getPixels(), height, width);

        for (const IntPixel_t &pix: maxRect) {
            const IntPoint &ll = pix.first;
            const IntPoint &ur = pix.second;

            Rectangle rect;
            rect.width_m = (ur.x_m - ll.x_m) * pixel_width;
            rect.height_m = (ur.y_m - ll.y_m) * pixel_height;

            double midX = 0.5 * (ur.x_m + ll.x_m);
            double midY = 0.5 * (ur.y_m + ll.y_m);
            rect.trafo_m = AffineTransformation(Vector_t(1, 0, (0.5 * width - midX) * pixel_width),
                                                Vector_t(0, 1, (midY - 0.5 * height) * pixel_height));

            pixmap->pixels_m.push_back(rect);
            pixmap->pixels_m.back().computeBoundingBox();
        }

        it += (arguments.getLengthConsumed() + 1);

        return true;
    }

    void Mask::print(int ident) {
        for (auto pix: pixels_m) pix.print(ident);
    }

    void Mask::apply(std::vector<std::shared_ptr<Base> > &bfuncs) {
        for (auto pix: pixels_m) pix.apply(bfuncs);
    }
}