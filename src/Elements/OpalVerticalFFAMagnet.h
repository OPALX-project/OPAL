//
// Header file for OpalVerticalFFAMagnet Element
//
// Copyright (c) 2019 Chris Rogers
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.

#ifndef OPAL_OPALVERTICALFFAMAGNET_H
#define OPAL_OPALVERTICALFFAMAGNET_H

#include "Elements/OpalElement.h"

/** OpalVerticalFFAMagnet provides the user interface for the VERTICALFFA object
 */
class OpalVerticalFFAMagnet : public OpalElement {
  public:
    /** enum maps string to integer value for UI definitions */
    enum {
        B0 = COMMON,
        FIELD_INDEX,
        WIDTH,
        MAX_HORIZONTAL_POWER,
        END_LENGTH,
        CENTRE_LENGTH,
        BB_LENGTH,
        HEIGHT_NEG_EXTENT,
        HEIGHT_POS_EXTENT,
        SIZE // size of the enum
    };

    /** Default constructor initialises UI parameters. */
    OpalVerticalFFAMagnet();

    /** Destructor does nothing */
    virtual ~OpalVerticalFFAMagnet();

    /** Inherited copy constructor */
    virtual OpalVerticalFFAMagnet *clone(const std::string &name);

    /** Fill in all registered attributes
     *
     *  Just calls fillRegisteredAttributes on the base class
     */
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /** Update the VerticalFFA with new parameters from UI parser */
    virtual void update();

  private:
    // Not implemented.
    OpalVerticalFFAMagnet(const OpalVerticalFFAMagnet &);
    void operator=(const OpalVerticalFFAMagnet &);

    // Clone constructor.
    OpalVerticalFFAMagnet(const std::string &name, OpalVerticalFFAMagnet *parent);

    static std::string docstring_m;
};

#endif // OPAL_OPALVERTICALFFAMAGNET_H

