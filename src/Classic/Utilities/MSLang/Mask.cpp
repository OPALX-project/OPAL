#include "Utilities/MSLang/Mask.h"
#include "Utilities/MSLang/ArgumentExtractor.h"
#include "Utilities/MSLang/matheval.h"
#include "Utilities/PortableBitmapReader.h"

#include <boost/regex.hpp>
#include <filesystem>

namespace {
    void updateCache(const std::vector<bool> &pixels, std::vector<int> &cache, int x, int width) {
        const int height = cache.size();
        for (int y = 0; y < height; ++ y) {
            if (pixels[x  + y * width]) {
                ++ cache[y];
            } else {
                cache[y] = 0;
            }
        }
    }
}

namespace mslang {
    long Mask::computeArea(const Mask::IntPoint &ll, const Mask::IntPoint &ur) const {
        if ((ur.x_m > ll.x_m) && (ur.y_m > ll.y_m))
            return (ur.x_m - ll.x_m) * (ur.y_m - ll.y_m);

        return 0;
    }

    std::pair<Mask::IntPoint, Mask::IntPoint> Mask::findMaximalRectangle(const std::vector<bool> &pixels,
                                                                         int M, /* height */
                                                                         int N /* width */) const {

        // This algorithm was presented in
        // http://www.drdobbs.com/database/the-maximal-rectangle-problem/184410529
        // by David Vandevoorde, April 01, 1998

        int bestArea = 0;
        IntPoint bestLL(0, 0), bestUR(0,0);
        std::vector<int> cache(M, 0);
        std::stack<std::pair<int, int> > stack;
        for (int x = N - 1; x >= 0; -- x) {
            updateCache(pixels, cache, x, N);
            int width = 0;
            for (int y = 0; y < M; ++ y) {
                if (cache[y] > width) {
                    stack.push(std::make_pair(y, width));
                    width = cache[y];
                } else if (cache[y] < width) {
                    int y0, w0;
                    do {
                        std::tie(y0, w0) = stack.top();
                        stack.pop();
                        if (width * (y - y0) > bestArea) {
                            bestLL.x_m = x; bestLL.y_m = y0;
                            bestUR.x_m = x + width; bestUR.y_m = y;
                            bestArea = width * (y - y0);
                        }
                        width = w0;
                    } while (!stack.empty() && cache[y] < width);
                    width = cache[y];
                    if (width != 0) {
                        stack.push(std::make_pair(y0, width));
                    }
                }
            }

            if (!stack.empty()) {
                int y0, w0;
                std::tie(y0, w0) = stack.top();
                stack.pop();
                if (width * (N - y0) > bestArea) {
                    bestLL.x_m = x; bestLL.y_m = y0;
                    bestUR.x_m = x + width; bestUR.y_m = N;
                    bestArea = width * (N - y0);
                }
            }
        }

        return std::make_pair(bestLL, bestUR);
    }

    std::vector<Mask::IntPixel_t> Mask::minimizeNumberOfRectangles(std::vector<bool> pixels,
                                                                   int M,/* height */
                                                                   int N /* width */)
    {
        std::vector<IntPixel_t> rectangles;

        long maxArea = 0;
        while (true) {
            IntPixel_t pix = findMaximalRectangle(pixels, M, N);
            long area = computeArea(pix.first, pix.second);
            if (area > maxArea) maxArea = area;
            if (area <= 1) {
                break;
            }

            rectangles.push_back(pix);

            for (int y = pix.first.y_m; y < pix.second.y_m; ++ y) {
                int idx = y * N + pix.first.x_m;
                for (int x = pix.first.x_m; x < pix.second.x_m; ++ x, ++ idx) {
                    pixels[idx] = false;
                }
            }
        }
        int idx = 0;
        for (int y = 0; y < M; ++ y) {
            for (int x = 0; x < N; ++ x, ++idx) {
                if (pixels[idx]) {
                    IntPoint ll(x, y);
                    IntPoint ur(x, y);
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

        if (!std::filesystem::exists(filename)) {
            ERRORMSG("file '" << filename << "' doesn't exists" << endl);
            return false;
        }

        PortableBitmapReader reader(filename);
        int width = reader.getWidth();
        int height = reader.getHeight();

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