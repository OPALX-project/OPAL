// ------------------------------------------------------------------------
// $RCSfile: TrackEnd.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: TrackEnd
//   The class for the OPAL ENDTRACK command.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/TrackEnd.h"
#include "Track/Track.h"
#include "Track/TrackParser.h"


// Class TrackEnd
// ------------------------------------------------------------------------

TrackEnd::TrackEnd():
    Action(0, "ENDTRACK",
           "The \"ENDTRACK\" sub-command stops tracking.")
{}


TrackEnd::TrackEnd(const std::string &name, TrackEnd *parent):
    Action(name, parent)
{}


TrackEnd::~TrackEnd()
{}


TrackEnd *TrackEnd::clone(const std::string &name) {
    return new TrackEnd(name, this);
}


void TrackEnd::execute() {
    Track::block->parser.stop();
}
