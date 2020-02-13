#include "Utilities/Util.h"
#include "Physics/Physics.h"
#include "OPALrevision.h"

#include <boost/regex.hpp>

#include <queue>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>

namespace Util {
    std::string getGitRevision() {
        return std::string(GIT_VERSION);
    }

#define erfinv_a3 -0.140543331
#define erfinv_a2 0.914624893
#define erfinv_a1 -1.645349621
#define erfinv_a0 0.886226899

#define erfinv_b4 0.012229801
#define erfinv_b3 -0.329097515
#define erfinv_b2 1.442710462
#define erfinv_b1 -2.118377725
#define erfinv_b0 1

#define erfinv_c3 1.641345311
#define erfinv_c2 3.429567803
#define erfinv_c1 -1.62490649
#define erfinv_c0 -1.970840454

#define erfinv_d2 1.637067800
#define erfinv_d1 3.543889200
#define erfinv_d0 1

    double erfinv (double x) // inverse error function
    {
        double x2, r, y;
        int  sign_x;

        if (x < -1 || x > 1)
            return NAN;

        if (x == 0)
            return 0;

        if (x > 0)
            sign_x = 1;
        else {
            sign_x = -1;
            x = -x;
        }

        if (x <= 0.7) {

            x2 = x * x;
            r =
                x * (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
            r /= (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 +
                  erfinv_b1) * x2 + erfinv_b0;
        }
        else {
            y = sqrt (-log ((1 - x) / 2));
            r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
            r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
        }

        r = r * sign_x;
        x = x * sign_x;

        r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
        r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));

        return r;
    }

#undef erfinv_a3
#undef erfinv_a2
#undef erfinv_a1
#undef erfinv_a0

#undef erfinv_b4
#undef erfinv_b3
#undef erfinv_b2
#undef erfinv_b1
#undef erfinv_b0

#undef erfinv_c3
#undef erfinv_c2
#undef erfinv_c1
#undef erfinv_c0

#undef erfinv_d2
#undef erfinv_d1
#undef erfinv_d0

        Vector_t getTaitBryantAngles(Quaternion rotation, const std::string &/*elementName*/) {
        Quaternion rotationBAK = rotation;

        // y axis
        Vector_t tmp = rotation.rotate(Vector_t(0, 0, 1));
        tmp(1) = 0.0;
        // tmp /= euclidean_norm(tmp);
        double theta = fmod(atan2(tmp(0), tmp(2)) + Physics::two_pi, Physics::two_pi);

        Quaternion rotTheta(cos(0.5 * theta), 0, sin(0.5 * theta), 0);
        rotation = rotTheta.conjugate() * rotation;

        // x axis
        tmp = rotation.rotate(Vector_t(0, 0, 1));
        tmp(0) = 0.0;
        tmp /= euclidean_norm(tmp);
        double phi = fmod(atan2(-tmp(1), tmp(2)) + Physics::two_pi, Physics::two_pi);

        Quaternion rotPhi(cos(0.5 * phi), sin(0.5 * phi), 0, 0);
        rotation = rotPhi.conjugate() * rotation;

        // z axis
        tmp = rotation.rotate(Vector_t(1, 0, 0));
        tmp(2) = 0.0;
        tmp /= euclidean_norm(tmp);
        double psi = fmod(atan2(tmp(1), tmp(0)) + Physics::two_pi, Physics::two_pi);

        return Vector_t(theta, phi, psi);
    }

    std::string toUpper(const std::string &str) {
        std::string output = str;
        std::transform(output.begin(), output.end(), output.begin(), [](const char c) { return toupper(c);});

        return output;
    }

    KahanAccumulation::KahanAccumulation():
        sum(0.0),
        correction(0.0)
    { }


    KahanAccumulation& KahanAccumulation::operator+=(double value) {
        long double y = value - this->correction;
        long double t = this->sum + y;
        this->correction = (t - this->sum) - y;
        this->sum = t;
        return *this;
    }

    /** \brief
     *  rewind the SDDS file such that the spos of the last step is less or equal to maxSPos
     */
    unsigned int rewindLinesSDDS(const std::string &fileName, double maxSPos, bool checkForTime) {
        if (Ippl::myNode() > 0) return 0;

        std::fstream fs(fileName.c_str(), std::fstream::in);
        if (!fs.is_open()) return 0;

        std::string line;
        std::queue<std::string> allLines;
        unsigned int numParameters = 0;
        unsigned int numColumns = 0;
        unsigned int sposColumnNr = 0;
        unsigned int timeColumnNr = 0;
        double spos, time = 0.0;
        double lastTime = -1.0;

        boost::regex parameters("&parameter");
        boost::regex column("&column");
        boost::regex data("&data");
        boost::regex end("&end");
        boost::regex name("name=([a-zA-Z0-9\\$_]+)");
        boost::smatch match;

        std::istringstream linestream;

        while (getline(fs, line)) {
            allLines.push(line);
        }
        fs.close();


        fs.open (fileName.c_str(), std::fstream::out);

        if (!fs.is_open()) return 0;

        do {
            line = allLines.front();
            allLines.pop();
            fs << line << "\n";
            if (boost::regex_search(line, match, parameters)) {
                ++numParameters;
                while (!boost::regex_search(line, match, end)) {
                    line = allLines.front();
                    allLines.pop();
                    fs << line << "\n";
                }
            } else if (boost::regex_search(line, match, column)) {
                ++numColumns;
                while (!boost::regex_search(line, match, name)) {
                    line = allLines.front();
                    allLines.pop();
                    fs << line << "\n";
                }
                if (match[1] == "s") {
                    sposColumnNr = numColumns;
                }
                if (match[1] == "t") {
                    timeColumnNr = numColumns;
                }
                while (!boost::regex_search(line, match, end)) {
                    line = allLines.front();
                    allLines.pop();
                    fs << line << "\n";
                }
            }
        } while (!boost::regex_search(line, match, data));

        while (!boost::regex_search(line, match, end)) {
            line = allLines.front();
            allLines.pop();
            fs << line << "\n";
        }

        for (unsigned int i = 0; i < numParameters; ++ i) {
            fs << allLines.front() << "\n";
            allLines.pop();
        }

        while (allLines.size() > 0) {
            line = allLines.front();

            linestream.str(line);
            if (checkForTime) {
                for (unsigned int i = 0; i < timeColumnNr; ++ i) {
                    linestream >> time;
                }
            }

            linestream.str(line);
            for (unsigned int i = 0; i < sposColumnNr; ++ i) {
                linestream >> spos;
            }

            if ((spos - maxSPos) > 1e-20 * Physics::c) break;

            allLines.pop();

            if (!checkForTime || (time - lastTime) > 1e-20)
                fs << line << "\n";

            lastTime = time;
        }

        fs.close();

        if (allLines.size() > 0)
            INFOMSG(level2 << "rewind " + fileName + " to " + std::to_string(maxSPos) << " m" << endl);

        return allLines.size();
    }

    /*
   base64.cpp and base64.h

   Copyright (C) 2004-2008 René Nyffenegger

   This source code is provided 'as-is', without any express or implied
   warranty. In no event will the author be held liable for any damages
   arising from the use of this software.

   Permission is granted to anyone to use this software for any purpose,
   including commercial applications, and to alter it and redistribute it
   freely, subject to the following restrictions:

   1. The origin of this source code must not be misrepresented; you must not
      claim that you wrote the original source code. If you use this source code
      in a product, an acknowledgment in the product documentation would be
      appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
      misrepresented as being the original source code.

   3. This notice may not be removed or altered from any source distribution.

   René Nyffenegger rene.nyffenegger@adp-gmbh.ch

    */

    static const std::string base64_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                            "abcdefghijklmnopqrstuvwxyz"
                                            "0123456789+/";


    static inline bool is_base64(unsigned char c) {
        return (isalnum(c) || (c == '+') || (c == '/'));
    }

    std::string base64_encode(const std::string &string_to_encode) {
        const char* bytes_to_encode = string_to_encode.c_str();
        unsigned int in_len = string_to_encode.size();
        std::string ret;
        int i = 0;
        int j = 0;
        unsigned char char_array_3[3];
        unsigned char char_array_4[4];

        while (in_len--) {
            char_array_3[i++] = *(bytes_to_encode++);
            if (i == 3) {
                char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
                char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
                char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
                char_array_4[3] = char_array_3[2] & 0x3f;

                for(i = 0; (i <4) ; i++)
                    ret += base64_chars[char_array_4[i]];
                i = 0;
            }
        }

        if (i)
            {
                for(j = i; j < 3; j++)
                    char_array_3[j] = '\0';

                char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
                char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
                char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
                char_array_4[3] = char_array_3[2] & 0x3f;

                for (j = 0; (j < i + 1); j++)
                    ret += base64_chars[char_array_4[j]];

                while((i++ < 3))
                    ret += '=';

            }

        return ret;
    }

    std::string base64_decode(std::string const& encoded_string) {
        int in_len = encoded_string.size();
        int i = 0;
        int j = 0;
        int in_ = 0;
        unsigned char char_array_4[4], char_array_3[3];
        std::string ret;

        while (in_len-- && ( encoded_string[in_] != '=') && is_base64(encoded_string[in_])) {
            char_array_4[i++] = encoded_string[in_]; in_++;
            if (i ==4) {
                for (i = 0; i <4; i++)
                    char_array_4[i] = base64_chars.find(char_array_4[i]);

                char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
                char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
                char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

                for (i = 0; (i < 3); i++)
                    ret += char_array_3[i];
                i = 0;
            }
        }

        if (i) {
            for (j = i; j <4; j++)
                char_array_4[j] = 0;

            for (j = 0; j <4; j++)
                char_array_4[j] = base64_chars.find(char_array_4[j]);

            char_array_3[0] = (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
            char_array_3[1] = ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
            char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

            for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
        }

        return ret;
    }
}
