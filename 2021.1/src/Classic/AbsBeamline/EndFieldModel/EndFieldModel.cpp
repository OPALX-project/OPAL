/*
 *  Copyright (c) 2017, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to
 *     endorse or promote products derived from this software without specific
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include <algorithm>

#include "AbsBeamline/EndFieldModel/EndFieldModel.h"

namespace endfieldmodel {

bool GreaterThan(std::vector<int> v1, std::vector<int> v2) {
  size_t n1(v1.size()), n2(v2.size());
  for (size_t i = 0; i < n1 && i < n2; ++i) {
    if (v1[n1-1-i] > v2[n2-1-i]) return true;
    if (v1[n1-1-i] < v2[n2-1-i]) return false;
  }
  return false;
}

std::vector< std::vector<int> > CompactVector(
                              std::vector< std::vector<int> > vec) {
  // first sort the list
  std::sort(vec.begin(), vec.end(), GreaterThan);
  // now look for n = n+1
  for (size_t j = 0; j < vec.size()-1; ++j) {
    while (j < vec.size()-1 && IterableEquality(
                              vec[j].begin()+1, vec[j].end(),
                              vec[j+1].begin()+1, vec[j+1].end()) ) {
      vec[j][0] += vec[j+1][0];
      vec.erase(vec.begin()+j+1);
    }
  }
  return vec;
}

}

