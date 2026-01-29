/*
 *  Copyright (c) 2014, Chris Rogers
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

#ifndef _CLASSIC_SRC_ALGORITHMS_ABSTRACTTIMEDEPENDENCE_H_
#define _CLASSIC_SRC_ALGORITHMS_ABSTRACTTIMEDEPENDENCE_H_

#include <string>
#include <map>
#include <memory>

/** @class AbstractTimeDependence
 * 
 *  Time dependence abstraction for field parameters that vary slowly with time;
 *  for example, in RF cavities for synchrotrons the RF frequency varies to
 *  match the time-of-flight of the particles
 * 
 *  This base class stores a mapping of string to time dependence. At Visit
 *  time, we do a map look-up to assign the appropriate TimeDependence to the
 *  relevant field Elements.
 */
class AbstractTimeDependence {
  // Unit tests for this class are in with PolynomialTimeDependence
  public:
    /** Destructor does nothing */
    virtual ~AbstractTimeDependence() {}
    /** Inheritable copy constructor
     *
     *  @returns new AbstractTimeDependence that is a copy of this. User owns
     *  returned memory.
     */
     virtual AbstractTimeDependence* clone() = 0;

    /** getValue(time) returns the value as a function of time.
     *
     *  This could represent RF voltage or frequency, magnetic field strength 
     *  etc
     */
    virtual double getValue(double time) = 0;

    /** Look up the time dependence that has a given name
     *  
     *  @param name name of the time dependence
     * 
     *  @returns shared_ptr to the appropriate time dependence.
     *  @throws GeneralClassicException if name is not recognised
     */
    static std::shared_ptr<AbstractTimeDependence>
                                            getTimeDependence(std::string name);

    /** Add a value to the lookup table
     *  
     *  @param name name of the time dependence. If name already exists in the
     *  map, it is overwritten with the new value.
     *  @param time_dep shared_ptr to the time dependence.
     */
    static void setTimeDependence(std::string name,
                              std::shared_ptr<AbstractTimeDependence> time_dep);

    /** Get the name corresponding to a given time_dep
     *
     *  @param time_dep time dependence to lookup
     *
     *  @returns name corresponding to the time dependence. Note that this
     *  just does a dumb loop over the stored map values; so O(N).
     *  @throws GeneralClassicException if time_dep is not recognised
     */
    static std::string getName
                             (std::shared_ptr<AbstractTimeDependence> time_dep);

  private:
    static std::map<std::string, std::shared_ptr<AbstractTimeDependence> > td_map;
};

#endif
