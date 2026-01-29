#ifndef OPAL_PlaceRep_HH
#define OPAL_PlaceRep_HH

// ------------------------------------------------------------------------
// $RCSfile: PlaceRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: PlaceRep
//
// ------------------------------------------------------------------------
//
// $Date: 2001/08/24 19:35:03 $
// $Author: jsberg $
//
// ------------------------------------------------------------------------

#include <string>
#include <vector>
#include <iosfwd>//tream>

class FlaggedElmPtr;


// Class PlaceRep
// ------------------------------------------------------------------------
/// Representation of a place within a beam line or sequence.
//  This representation holds data to keep track wether the place was
//  already passed when tracking through the line.  The place may be
//  one of the following:
//  [ul]
//  [li]"#S", the start of the beam line.
//  [li]"#E", the end of the beam line.
//  [li]"SELECTED", all selected positions (flag set in FlaggedElmPtr).
//  [li]name[occurrence]::name[occurrence]..., denoting a specific place.
//  [/ul]

class PlaceRep {

public:

    /// Default constructor.
    //  Constructs undefined place.
    PlaceRep();

    /// Construct default place.
    //  Used for places like "#S", "#E", "SELECTED".
    explicit PlaceRep(const std::string &);

    PlaceRep(const PlaceRep &);
    ~PlaceRep();
    const PlaceRep &operator=(const PlaceRep &);

    /// Add a name/occurrence pair.
    void append(const std::string &, int occur);

    /// Initialise data for search.
    //  Sets the internal state to the beginning of the line.
    void initialize();

    /// Enter an element or line.
    //  Sets the internal state to active, when the specified place is
    //  entered.
    void enter(const FlaggedElmPtr &) const;

    /// Leave an element or line.
    //  Sets the internal state to inactive, when the specified place is
    //  left.
    void leave(const FlaggedElmPtr &) const;

    /// Return status.
    //  Returns true, if we are at the specified place.
    bool isActive() const;

    /// Return select flag.
    //  Returns true, if the place has the value "SELECTED".
    bool isSelected() const;

    /// Print in input format.
    void print(std::ostream &os) const;

private:

    typedef std::vector<std::pair<std::string, int> > Data;
    Data data;

    // If true, this place refers to "SELECTED"
    bool is_selected;

    // If true, we are at this position.
    mutable bool status;

    // Ancillary variable, keeps track of levels seen.
    mutable Data::size_type seen;
};


inline std::ostream &operator<<(std::ostream &os, const PlaceRep &data) {
    data.print(os);
    return os;
}

#endif // OPAL_PlaceRep_HH
