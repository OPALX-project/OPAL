#include <map>
#include <string>

#include "changes.h"

namespace Versions {
    std::map<unsigned int, std::string> changes;

    void fillChanges() {
        if (changes.size() > 0) return;

        changes.insert({105,
                    "* The normalization of the 2-dimensional field maps has changed.\n"
                    "  Instead of normalizing with the overall maximum value of longitudinal\n"
                    "  component Opal now uses the maximum value on axis.\n"
                    "\n"
                    "* The parser has been modified to check the type of all variables. All real \n"
                    "  variables have to be prefixed with the keyword REAL.\n"
                    });

        changes.insert({109,
                    "* The attribute BFREQ of the command BEAM is now in MHz instead of Hz\n"
                    "\n"
                    "* OPAL-T: Beamlines are now placed in 3-dimensional space. Make sure that\n"
                    "  you use apertures to limit the range of the elements. Default aperture \n"
                    "  has circular shape with diameter 1 meter.\n"
                    "\n"
                    "* OPAL-T: Beamlines containing a cathode have to have a SOURCE element to\n"
                    "  indicate this fact\n"
                    "\n"
                    "* OPAL-T: The design energy of dipoles is now expected to be in MeV instead\n"
                    "  of eV.\n"
                    "\n"
                    "* OPAL-T: The attribute 'ROTATION' of RBEND and SBEND has been replaced\n"
                    "  by 'PSI'. Can be applyied to all elements to rotate them.\n"
                    "\n"
                    "* The attribute DISTRIBUTION of the command DISTRIBUTION has been renamed to\n"
                    "  TYPE.\n"
                    "\n"
                    "* The meaning of OFFSETZ of the command DISTRIBUTION has changed. It now \n"
                    "  indicates a shift of the particle bunch relative to the reference particle.\n"
                    "  Use the ZSTART attribute of the TRACK command to start the simulation at a \n"
                    "  position z > 0.\n"
                    "\n"
                    "* The string indicating the orientation (sofar always XYZ) of 3D fieldmaps has\n"
                    "  been dropped.\n"
                    });
    }
}