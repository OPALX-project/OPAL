//
// Struct SumErrSqRadialPeak
//   A simple expression computing the sum of all peak errors (given as
//   first and second argument) for a range of peaks (third argument and fourth argument)
//   according to
//
//   \f[
//     result = \frac{1}{n} * \sqrt{\sum_{i=start}^end (measurement_i - value_i)^2}
//   \f]
//
// Copyright (c) 2017 - 2018, Matthias Frey, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// Implemented as part of the PhD thesis
// "Precise Simulations of Multibunches in High Intensity Cyclotrons"
// and the paper
// "Matching of turn pattern measurements for cyclotrons using multiobjective optimization"
// (https://doi.org/10.1103/PhysRevAccelBeams.22.064602)
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
#include "Expression/SumErrSqRadialPeak.h"

const std::string name = "SumErrSqRadialPeak";
