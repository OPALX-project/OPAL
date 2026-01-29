/*

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

 ----------------------------------------------------------------------

*/

#include <string>
#include <vector>

// To do:
// - can we sort less often or reduce/optimise dominance checks?
// - should we use FPL's data structure?
// - two changes in read.c
// - heuristics

// opt:  0 = basic, 1 = sorting, 2 = slicing to 2D, 3 = slicing to 3D

//For some reason opt=3 doesnt work for some of the test data
//So keep opt = 2
#define hyper_opt 2

namespace Hypervolume {

    double FromFile(std::string file, const std::vector<double>& referencePoint);
    //Room for more functions to compute volumes
    //without accessing external files
}
