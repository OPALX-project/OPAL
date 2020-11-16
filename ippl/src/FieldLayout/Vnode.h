y//
// Vnodes really have very little information.
// It knows its domain and what processor it resides on.
//
// Also, it has a global integer index for the vnode (useful with more recent
// FieldLayouts which store a logical "array" of vnodes; user specifies numbers
// of vnodes along each direction). Classes or user codes that use Vnode are
// responsible for setting and managing the values of this index; if unset, it
// has the value -1.
//
// Copyright (c) 2003 - 2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
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

#ifndef VNODE_H
#define VNODE_H

#include "Index/NDIndex.h"
#include "Utility/RefCounted.h"

#include <iostream>

// forward declarations
template <unsigned Dim>
class Vnode;
template <unsigned Dim>
std::ostream& operator<<(std::ostream&, const Vnode<Dim>&);

//----------------------------------------------------------------------

template <unsigned Dim>
class Vnode : public RefCounted {
private:
    NDIndex<Dim> Domain;
    int Node;
    int vnode_m;  // Global vnode ID number (between 0 and nvnodes - 1)

public:
    // Null ctor does nothing.
    Vnode() {
    }

    // Normal ctor:
    Vnode(const NDIndex<Dim>& domain, int node, int vnode = -1)
        : Domain(domain), Node(node), vnode_m(vnode) {
    }

    // Copy ctor:
    Vnode(const Vnode<Dim>& vn) : Domain(vn.Domain), Node(vn.Node), vnode_m(vn.vnode_m) {
    }

    // operator= to copy one vnode into another
    Vnode<Dim>& operator=(const Vnode<Dim>& vn) {
        Domain  = vn.Domain;
        Node    = vn.Node;
        vnode_m = vn.vnode_m;
        return *this;
    }

    int getNode() const {
        return Node;
    }
    int getVnode() const {
        return vnode_m;
    }
    const NDIndex<Dim>& getDomain() const {
        return Domain;
    }

    // put data into a message to send to another node
    Message& putMessage(Message& m) const {
        Domain.putMessage(m);
        m.put(Node);
        m.put(vnode_m);
        return m;
    }

    // get data out from a message
    Message& getMessage(Message& m) {
        Domain.getMessage(m);
        m.get(Node);
        m.get(vnode_m);
        return m;
    }
};

template <unsigned Dim>
inline std::ostream& operator<<(std::ostream& out, const Vnode<Dim>& v) {
    out << "Node = " << v.getNode() << " ; vnode_m = " << v.getVnode()
        << " ; Domain = " << v.getDomain();
    return out;
}

#endif
