#ifndef OPAL_NILTRACKER_H
#define OPAL_NILTRACKER_H

//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#define NIL_VISITELEMENT(elem) virtual void visit##elem(const elem &) { }

#include "Algorithms/Tracker.h"

class BMultipoleField;
template <class T, unsigned Dim>
class PartBunchBase;
class AlignWrapper;
class BeamBeam;
class BeamStripping;
class CCollimator;
class Corrector;
class CyclotronValley;
class Degrader;
class Diagnostic;
class Drift;
class ElementBase;
class FlexibleCollimator;
class Lambertson;
class Marker;
class Monitor;
class Multipole;
class ParallelPlate;
class Probe;
class RBend;
class RFCavity;
class RFQuadrupole;
class SBend;
class Separator;
class Septum;
class Solenoid;
class TravelingWave;

class NilTracker: public Tracker {

public:
    /// Constructor.
    explicit NilTracker(const Beamline &beamline,
                        const PartData &reference,
                        bool revBeam,
                        bool revTrack);

    virtual ~NilTracker();

    NIL_VISITELEMENT(AlignWrapper)
    NIL_VISITELEMENT(Beamline)
    NIL_VISITELEMENT(BeamBeam)
    NIL_VISITELEMENT(BeamStripping)
    NIL_VISITELEMENT(CCollimator)
    NIL_VISITELEMENT(Corrector)
    NIL_VISITELEMENT(CyclotronValley)
    NIL_VISITELEMENT(Degrader)
    NIL_VISITELEMENT(Diagnostic)
    NIL_VISITELEMENT(Drift)
    NIL_VISITELEMENT(FlexibleCollimator)
    NIL_VISITELEMENT(Lambertson)
    NIL_VISITELEMENT(Marker)
    NIL_VISITELEMENT(Monitor)
    NIL_VISITELEMENT(Multipole)
    NIL_VISITELEMENT(ParallelPlate)
    NIL_VISITELEMENT(Probe)
    NIL_VISITELEMENT(RBend)
    NIL_VISITELEMENT(RFCavity)
    NIL_VISITELEMENT(RFQuadrupole)
    NIL_VISITELEMENT(SBend)
    NIL_VISITELEMENT(Separator)
    NIL_VISITELEMENT(Septum)
    NIL_VISITELEMENT(Solenoid)
    NIL_VISITELEMENT(TravelingWave)

    virtual void execute();

private:

    NilTracker();
    NilTracker(const NilTracker &);

    void operator=(const NilTracker &);
};

#endif // OPAL_NILTRACKER_H