#ifndef OPAL_Directory_HH
#define OPAL_Directory_HH

// ------------------------------------------------------------------------
// $RCSfile: Directory.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Directory
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:34 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "MemoryManagement/Pointer.h"
#include <functional>
#include <map>
#include <string>

class Object;


// Class Directory
// ------------------------------------------------------------------------

typedef std::map<std::string, Pointer<Object>, std::less<std::string> >
ObjectDir;

/// A map of string versus pointer to Object.
//  Used as the directory for OPAL objects.

class Directory {

public:

    /// Constructor.
    //  Build empty directory.
    Directory();

    ~Directory();

    /// First object in alphabetic order of name
    //  (Version for non-constant directory).
    ObjectDir::iterator begin();

    /// First object in alphabetic order of name
    //  (Version for constant directory).
    ObjectDir::const_iterator begin() const;

    /// Last object in alphabetic order of name
    //  (Version for non-constant directory).
    ObjectDir::iterator end();

    /// Last object in alphabetic order of name
    //  (Version for constant directory).
    ObjectDir::const_iterator end() const;

    /// Delete all entries.
    void erase();

    /// Remove existing entry.
    //  The entry is identified by [b]name[/b].
    void erase(const std::string &name);

    /// Find entry.
    //  The entry is identified by [b]name[/b].
    //  If the entry [b]name[/b] does not exist, return [b]NULL[/b].
    Object *find(const std::string &name) const;

    /// Define new object.
    //  Insert new object in directory.
    //  If the entry [b]name[/b] exists already, it is removed.
    void insert(const std::string &name, Object *newObject);

private:

    // The directory of objects.
    ObjectDir dir;
};

#endif // OPAL_Directory_HH
