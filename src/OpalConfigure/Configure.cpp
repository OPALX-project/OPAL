// ------------------------------------------------------------------------
// $RCSfile: Configure.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Namespace: Configure
//   Contains methods for configuring the OPAL-9 program.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/05/03 12:40:49 $
// $Author: opal $
//
// JMJ & JP adding Aperture and Split 18/4/2000
// ------------------------------------------------------------------------

#include "OpalConfigure/Configure.h"
#include "AbstractObjects/OpalData.h"

#include "Distribution/Distribution.h"

// Basic action commands.
#include "BasicActions/Call.h"
#include "BasicActions/Dump.h"
#include "BasicActions/DumpFields.h"
#include "BasicActions/DumpEMFields.h"
#include "BasicActions/Echo.h"
#include "BasicActions/Help.h"
#include "BasicActions/Option.h"
#include "BasicActions/Save.h"
#include "BasicActions/Select.h"
#include "BasicActions/Show.h"
#include "BasicActions/Stop.h"
#include "BasicActions/Quit.h"
#include "BasicActions/What.h"
#include "BasicActions/System.h"
#include "BasicActions/PSystem.h"
#include "BasicActions/Title.h"
#include "BasicActions/Value.h"

// Macro command.
#include "OpalParser/MacroCmd.h"

// Physics action commands.
#include "PhysicsActions/Dynamic.h"
#include "PhysicsActions/MakeSequence.h"
#include "PhysicsActions/SetIntegrator.h"
#include "PhysicsActions/Static.h"

// Commands introducing a special mode.
#include "Editor/EditCmd.h"
#include "Match/MatchCmd.h"
#include "Track/TrackCmd.h"

// Table-related commands.
#include "Structure/Beam.h"
#include "Structure/FieldSolver.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/OpalWake.h"
#include "Structure/ParticleMatterInteraction.h"
#include "Utilities/OpalFilter.h"
#include "TrimCoils/OpalTrimCoil.h"
#include "Tables/AttList.h"
#include "Tables/Eigen.h"
#include "Tables/Envelope.h"
#include "Tables/Insertion.h"
#include "Tables/List.h"
#include "Tables/MatrixCmd.h"
#include "Tables/Micado.h"
#include "Tables/Period.h"
#include "Tables/Survey.h"
#include "Tables/ThreadAll.h"
#include "Tables/ThreadBpm.h"
#include "Tables/Twiss3.h"
#include "Aperture/Aperture.h"
#include "Aperture/Split.h"


// Value definitions commands.
#include "ValueDefinitions/BoolConstant.h"
#include "ValueDefinitions/RealConstant.h"
#include "ValueDefinitions/RealVariable.h"
#include "ValueDefinitions/RealVector.h"
#include "ValueDefinitions/StringConstant.h"

// Element commands.
#include "Elements/OpalBeamStripping.h"
#include "Elements/OpalCavity.h"
#include "Elements/OpalCCollimator.h"
#include "Elements/OpalCyclotron.h"
#include "Elements/OpalDrift.h"
#include "Elements/OpalECollimator.h"
#include "Elements/OpalFlexibleCollimator.h"
#include "Elements/OpalDegrader.h"
#include "Elements/OpalHKicker.h"
#include "Elements/OpalHMonitor.h"
#include "Elements/OpalInstrument.h"
#include "Elements/OpalKicker.h"
#include "Elements/OpalMarker.h"
#include "Elements/OpalMonitor.h"
#include "Elements/OpalMultipole.h"
#include "Elements/OpalMultipoleT.h"
#include "Elements/OpalMultipoleTStraight.h"
#include "Elements/OpalMultipoleTCurvedConstRadius.h"
#include "Elements/OpalMultipoleTCurvedVarRadius.h"
#include "Elements/OpalOctupole.h"
#include "Elements/OpalOffset/OpalLocalCartesianOffset.h"
#include "Elements/OpalOffset/OpalLocalCylindricalOffset.h"
#include "Elements/OpalOffset/OpalGlobalCartesianOffset.h"
#include "Elements/OpalOffset/OpalGlobalCylindricalOffset.h"
#include "Elements/OpalPepperPot.h"
#include "Elements/OpalPatch.h"
#include "Elements/OpalProbe.h"
#include "Elements/OpalQuadrupole.h"
#include "Elements/OpalPolynomialTimeDependence.h"
#include "Elements/OpalRBend.h"
#include "Elements/OpalRBend3D.h"
#include "Elements/OpalRCollimator.h"
#include "Elements/OpalSBend.h"
#include "Elements/OpalSBend3D.h"
#include "Elements/OpalScalingFFAMagnet.h"
#include "Elements/OpalSeparator.h"
#include "Elements/OpalSeptum.h"
#include "Elements/OpalSextupole.h"
#include "Elements/OpalSlit.h"
#include "Elements/OpalSolenoid.h"
#include "Elements/OpalSource.h"
#include "Elements/OpalSRot.h"
#include "Elements/OpalTravelingWave.h"
#include "Elements/OpalVKicker.h"
#include "Elements/OpalVMonitor.h"
//#include "Elements/OpalWire.h"
#include "Elements/OpalYRot.h"
#include "Elements/OpalParallelPlate.h"
#include "Elements/OpalCyclotronValley.h"
#include "Elements/OpalStripper.h"
#include "Elements/OpalRingDefinition.h"
#include "Elements/OpalVariableRFCavity.h"
#include "Elements/OpalVariableRFCavityFringeField.h"
#include "Elements/OpalVerticalFFAMagnet.h"

// Structure-related commands.
#include "Lines/Line.h"
#include "Lines/Sequence.h"

// Optimize command
#include "Optimize/OptimizeCmd.h"
#include "Optimize/DVar.h"
#include "Optimize/Objective.h"
#include "Optimize/Constraint.h"

// Sample command
#include "Sample/SampleCmd.h"
#include "Sample/OpalSample.h"

#include "changes.h"

// Modify these methods to add new commands.
// ------------------------------------------------------------------------

namespace {

    void makeActions() {
        OpalData *opal = OpalData::getInstance();
        opal->create(new Call());
        opal->create(new Dump());
        opal->create(new DumpFields());
        opal->create(new DumpEMFields());
        opal->create(new Echo());
        opal->create(new Dynamic());
        opal->create(new Eigen());
        opal->create(new Envelope());
        opal->create(new Help());
        opal->create(new EditCmd());
        opal->create(new List());
        opal->create(new MakeSequence());
        opal->create(new MatchCmd());
        opal->create(new MatrixCmd());
        opal->create(new Micado());
        opal->create(new Option());
        opal->create(new OptimizeCmd());
        opal->create(new SampleCmd());
        opal->create(new Save());
        opal->create(new Select());
        opal->create(new Show());
        opal->create(new SetIntegrator());
        opal->create(new Static());
        opal->create(new Stop());
        opal->create(new Quit());
        opal->create(new PSystem());
        opal->create(new System());
        opal->create(new ThreadAll());
        opal->create(new ThreadBpm());
        opal->create(new Title());
        opal->create(new TrackCmd());
        opal->create(new Twiss3());
        opal->create(new Aperture());
        opal->create(new MSplit());
        opal->create(new Value());
        opal->create(new What());
    }


    void makeDefinitions() {
        OpalData *opal = OpalData::getInstance();
        // Must create the value definitions first.
        opal->create(new BoolConstant());
        opal->create(new RealConstant());
        opal->create(new RealVariable());
        opal->create(new RealVector());
        opal->create(new StringConstant());

        opal->create(new AttList());
        opal->create(new Beam());
        opal->create(new FieldSolver());
        opal->create(new BoundaryGeometry());
        opal->create(new OpalWake());
        opal->create(new ParticleMatterInteraction());

        opal->create(new OpalFilter());
        opal->create(new OpalTrimCoil());

        opal->create(new Distribution());

        opal->create(new MacroCmd());
        opal->create(new Period());
        opal->create(new Insertion());
        opal->create(new Survey());

        opal->create(new DVar());
        opal->create(new Objective());
        opal->create(new Constraint());

        opal->create(new OpalSample());
    }


    void makeElements() {
        OpalData *opal = OpalData::getInstance();
        opal->create(new OpalBeamStripping());
        opal->create(new OpalCavity());
        opal->create(new OpalCCollimator());
        opal->create(new OpalCyclotron());
        opal->create(new OpalDrift());
        opal->create(new OpalECollimator());
        opal->create(new OpalFlexibleCollimator());
        opal->create(new OpalDegrader());
        opal->create(new OpalHKicker());
        opal->create(new OpalHMonitor());
        opal->create(new OpalInstrument());
        opal->create(new OpalKicker());
        opal->create(new OpalMarker());
        opal->create(new OpalMonitor());
        opal->create(new OpalMultipole());
        opal->create(new OpalMultipoleT());
        opal->create(new OpalMultipoleTStraight());
        opal->create(new OpalMultipoleTCurvedConstRadius());
        opal->create(new OpalMultipoleTCurvedVarRadius());
        opal->create(new OpalOctupole());
        opal->create(new OpalOffset::OpalLocalCartesianOffset());
//        opal->create(new OpalOffset::OpalLocalCylindricalOffset());
//        opal->create(new OpalOffset::OpalGlobalCartesianOffset());
//        opal->create(new OpalOffset::OpalGlobalCylindricalOffset());
        opal->create(new OpalPatch());
        opal->create(new OpalProbe());
        opal->create(new OpalPepperPot());
        opal->create(new OpalPolynomialTimeDependence());
        opal->create(new OpalQuadrupole());
        opal->create(new OpalRBend());
        opal->create(new OpalRBend3D());
        opal->create(new OpalRCollimator());
        opal->create(new OpalSBend());
        opal->create(new OpalSBend3D());
        opal->create(new OpalScalingFFAMagnet());
        opal->create(new OpalSeparator());
        opal->create(new OpalSeptum());
        opal->create(new OpalSextupole());
        opal->create(new OpalSlit());
        opal->create(new OpalSolenoid());
        opal->create(new OpalSource());
        opal->create(new OpalSRot());
        opal->create(new OpalTravelingWave());
        opal->create(new OpalVariableRFCavity());
        opal->create(new OpalVariableRFCavityFringeField());
        opal->create(new OpalVerticalFFAMagnet());
        opal->create(new OpalVKicker());
        opal->create(new OpalVMonitor());
        // opal->create(new OpalWire());
        opal->create(new OpalYRot());
        opal->create(new OpalParallelPlate());
        opal->create(new OpalCyclotronValley());
        opal->create(new OpalStripper());
        opal->create(new Line());
        opal->create(new Sequence());
        opal->create(new OpalRingDefinition());
    }
};

namespace Configure {
    void configure() {
        makeDefinitions();
        makeElements();
        makeActions();
        Versions::fillChanges();
    }
};
