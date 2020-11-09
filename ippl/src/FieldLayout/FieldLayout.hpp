//
// FielLayout template class
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

#include "FieldLayout/FieldLayout.h"
#include "Message/Communicate.h"
#include "Message/Message.h"
#include "Utility/IpplInfo.h"
#include "Utility/IpplStats.h"
#include "Utility/PAssert.h"

#include <cstdlib>
#include <limits>

//////////////////////////////////////////////////////////////////////
// Default constructor, which should only be used if you are going to
// call 'initialize' soon after (before using in any context)
template<unsigned Dim>
FieldLayout<Dim>::FieldLayout() {

    // just initialize basic things
    vnodesPerDirection_m = 0;

    // for this kind of construction, we just take it that the user is
    // requesting all parallel axes
    for (unsigned int dl=0; dl < Dim; ++dl) RequestedLayout[dl] = PARALLEL;
}

//////////////////////////////////////////////////////////////////////
// Destructor: Everything deletes itself automatically ... the base
// class destructors inform all the FieldLayoutUser's we're going away.
template<unsigned Dim>
FieldLayout<Dim>::~FieldLayout() {
    if (vnodesPerDirection_m != 0)
        delete [] vnodesPerDirection_m;
}

//////////////////////////////////////////////////////////////////////

// Initialization functions, only to be called by the user of FieldLayout
// objects when the FieldLayout was created using the default constructor;
// otherwise these are only called internally by the various non-default
// FieldLayout constructors:

//-----------------------------------------------------------------------------
// These specify only a total number of vnodes, allowing the constructor
// complete control on how to do the vnode partitioning of the index space:

template<unsigned Dim>
void
FieldLayout<Dim>::initialize(const NDIndex<Dim>& domain,
			     e_dim_tag *p, int vnodes) {
    setup(domain, p, vnodes);
}

//////////////////////////////////////////////////////////////////////

//
// Now, the meat of the construction.
//

//-----------------------------------------------------------------------------
//This setup() specifies only a total number of vnodes, taking complete control
//on how to do the vnode partitioning of the index space:

template<unsigned Dim>
void
FieldLayout<Dim>::setup(const NDIndex<Dim>& domain,
			e_dim_tag *userflags, int vnodes)
{
    // Find the number processors.
    int nprocs = Ippl::getNodes();
    int myproc = Ippl::myNode();

    // If the user didn't specify the number of vnodes, make it equal nprocs
    if (vnodes <= 0) vnodes = nprocs;

    Inform dbgmsg("FieldLayout::setup", INFORM_ALL_NODES);
    // dbgmsg << "*** Domain=" << domain << ", nprocs=" << nprocs;
    // dbgmsg << ", myproc=" << myproc << ", vnodes=" << vnodes << endl;

    // If the user did not specify parallel/serial flags then make all parallel.
    int parallel_count = 0;
    unsigned int flagdim = 0;
    long totparelems = 1;
    for (flagdim=0; flagdim < Dim; ++flagdim) {
        if (userflags == 0)
            RequestedLayout[flagdim] = PARALLEL;
        else
            RequestedLayout[flagdim] = userflags[flagdim];
        if (RequestedLayout[flagdim] == PARALLEL) {
            parallel_count++;
            totparelems *= domain[flagdim].length();
        }
    }
    // Make sure at least one of the parallel/serial flags is parallel
    PInsist(parallel_count>0,"At least one dimension of a FieldLayout must be PARALLEL!");

    // Check to see if we have too few elements to partition.  If so, reduced
    // the number of vnodes (if necessary) to just the number of elements along
    // parallel dims.
    if (totparelems < vnodes) {
        //dbgmsg << "Total parallel lengths = " << totparelems << "; reducing ";
        //dbgmsg << "vnodes from " << vnodes << " to " << totparelems << endl;
        vnodes = totparelems;
    }

    e_dim_tag *flags = RequestedLayout;

    // Recursively split the domain until we have generated all the domains.
    Domain = domain;
    NDIndex<Dim> *domains_c = new NDIndex<Dim>[vnodes];
    NDIndex<Dim> *copy_c    = new NDIndex<Dim>[vnodes];
    NDIndex<Dim> leftDomain ;

    // Start with the whole domain.
    domains_c[0] = domain;
    int v;
    unsigned int d=0;

    int v1,v2,rm,vtot,vl,vr;
    double a,lmax,len;
    vnodesPerDirection_m = 0;
    for (v=vnodes,rm=0;v>1;v/=2) { rm += (v % 2); }
    if (rm == 0) {

        // vnodes is a power of 2

        // For power-of-two vnodes, allocate storage for vnodesPerDirection_m:
        //don't--tjw    vnodesPerDirection_m = new unsigned[Dim];
        // Fill in 1 for starters; remains 1 for serial directions:
        //don't--tjw    for (int d2=0; d2<Dim; d2++) vnodesPerDirection_m[d2] = 1;

        for (v=1; v<vnodes; v*=2) {
            // Go to the next parallel dimension.
            while(flags[d] != PARALLEL) if(++d == Dim) d = 0;

            // Split all the current vnodes.
            int i,j;
            for (i=0, j=0; i<v; ++i, j+=2)
                // Split to the left and to the right, saving both.
                domains_c[i].split( copy_c[j] , copy_c[j+1] , d );
            // Copy back.
            std::copy(copy_c,copy_c+v*2, domains_c);

            //don't--tjw      vnodesPerDirection_m[d] *= 2; // Update to reflect split

            // On to the next dimension.
            if (++d == Dim) d = 0;
        }

    } else {

        vtot = 1; // count the number of vnodes to make sure that it worked
        // vnodes is not a power of 2 so we need to do some fancy splitting
        // sorry... this would be much cleaner with recursion
        /*
          The way this works is to recursively split on the longest dimension.
          Suppose you request 11 vnodes.  It will split the longest dimension
          in the ratio 5:6 and put the new domains in node 0 and node 5.  Then
          it splits the longest dimension of the 0 domain and puts the results
          in node 0 and node 2 and then splits the longest dimension of node 5
          and puts the results in node 5 and node 8. etc.
          The logic is kind of bizarre, but it works.
        */
        for (v=1; v<2*vnodes; ++v) {
            // kind of reverse the bits of v
            for (v2=v,v1=1;v2>1;v2/=2) { v1 = 2*v1+(v2%2); }
            vl = 0; vr = vnodes;
            while (v1>1) {
                if ((v1%2)==1) {
                    vl=vl+(vr-vl)/2;
                } else {
                    vr=vl+(vr-vl)/2;
                }
                v1/=2;
            }
            v2=vl+(vr-vl)/2;

            if (v2>vl) {
                a = v2-vl;
                a /= vr-vl;
                vr=v2;
                leftDomain=domains_c[vl];
                lmax=0;
                d=std::numeric_limits<unsigned int>::max();
                for (unsigned int dd=0;dd<Dim;++dd) {
                    if ( flags[dd] == PARALLEL ) {
                        if ((len = leftDomain[dd].length()) > lmax) {
                            lmax = len;
                            d = dd;
                        }
                    }
                }
                domains_c[vl].split( domains_c[vl] , domains_c[vr] , d , a);
                ++vtot;
            }
        }
        v=vtot;
    }
    // Make sure we had a power of two number of vnodes.
    PAssert_EQ( v, vnodes );

    // Now make the vnodes, using the domains just generated.
    // Some of them we store in the local list, others in the remote.
    ac_domain_vnodes *remote_ac = new ac_domain_vnodes( domain );
    typedef typename ac_gc_domain_vnodes::value_type vtype;
    Remotes_ac.insert( vtype(gc0(),remote_ac) );
    for (v=0; v<vnodes; ++v)
        {
            int p = (v*nprocs)/vnodes;
            bool nosplit = (domains_c[v].size() < 2);

            // Add v arg to Vnode constructor 3/19/98 --tjw:
            Vnode<Dim> *vnode = new Vnode<Dim>(domains_c[v], p, v);
            typedef typename ac_id_vnodes::value_type v1;
            typedef typename ac_domain_vnodes::value_type v2;
            if ( p==myproc )
                Local_ac.insert(v1(Unique::get(),vnode));
            else
                remote_ac->insert(v2(domains_c[v], vnode), nosplit);
        }

    // Delete the memory we allocated.
    delete [] domains_c;
    delete [] copy_c;

    // Calculate the widths of all the vnodes
    calcWidths();
}

//-----------------------------------------------------------------------------
// Return number of vnodes along a direction.
template<unsigned Dim>
unsigned FieldLayout<Dim>::
getVnodesPerDirection(unsigned dir) {

    // If not set, then not appropriate to be calling this function. Example:
    // now-power-of-two-constructed FieldLayout which did *not* specify numbers
    // of vnodes along each direction:
    PAssert(vnodesPerDirection_m);

    return(vnodesPerDirection_m[dir]);
}

//----------------------------------------------------------------------
template<unsigned Dim>
FieldLayout<Dim>::FieldLayout(const NDIndex<Dim>& domain,
			      const NDIndex<Dim>* idx_begin,
			      const NDIndex<Dim>* idx_end)
    : Domain(domain)
{

    // we have one more FieldLayout, indicate this
    //INCIPPLSTAT(incFieldLayouts);

    // for this kind of construction, we just take it that the user is
    // requesting all parallel axes
    for (unsigned int dl=0; dl < Dim; ++dl) RequestedLayout[dl] = PARALLEL;

    // Build Vnodes for each of the local domains.
    vnodesPerDirection_m = 0;
    int mynode = Ippl::Comm->myNode();
    for (const NDIndex<Dim>* p=idx_begin; p!=idx_end; ++p)
        {
            typedef typename ac_id_vnodes::value_type v1;
            Local_ac.insert(v1(Unique::get(), new Vnode<Dim>(*p,mynode)));
        }

    // Everybody broadcasts their new local domains to everybody.
    // Build a message with the local domains and the id's for them.
    Message* bcast_mess = new Message;
    int count = idx_end - idx_begin;
    bcast_mess->put(count);
    // ::putMessage(*bcast_mess,idx_begin,idx_end);
    const NDIndex<Dim>* indxmsg = idx_begin;
    while (indxmsg != idx_end)
        (indxmsg++)->putMessage(*bcast_mess);
    // Send it to everybody except yourself, and record the number sent.
    int tag = Ippl::Comm->next_tag( F_REPARTITION_BCAST_TAG , F_TAG_CYCLE );
    int node_count = Ippl::Comm->broadcast_others(bcast_mess,tag);
    // Create the container for the remote vnodes.
    ac_domain_vnodes *remote_ac = new ac_domain_vnodes( Domain );
    typedef typename ac_gc_domain_vnodes::value_type v1;
    Remotes_ac.insert(v1(gc0(),remote_ac) );

    // Loop until we receive a message from each other node.
    while ((--node_count)>=0)
        {
            // Receive a broadcast message from any node.
            int other_node = COMM_ANY_NODE;
            Message *recv_mess = Ippl::Comm->receive_block(other_node,tag);
            PAssert(recv_mess);
            // Extract the number of vnodes coming in.
            int count = 0;
            recv_mess->get(count);
            // Now get the domains for the vnodes from the message.
            NDIndex<Dim> p;
            while (count-- > 0) {
                p.getMessage(*recv_mess);
                Vnode<Dim> *vnode = new Vnode<Dim>(p, other_node);
                // For domains of extent 1 in all directions, must call
                // DomainMap::insert() with the noSplit flag set to true, to avoid an
                // assertion failure in Index::split. Also do this for domains of
                // extent ZERO. These are not supposed to happen in IPPL, but when the
                // do when something like BinaryRepartition creates them, accomodating
                // them here will prevent subsequent code hangs or other errors and
                // allow recovery such as rejection of the result of BinaryRepartition
                // (discard the FieldLayout with zero-size domains and go on with the
                // original layout unchanged). (tjw)
                //	if (p.size() == 1) {
                if (p.size() <= 1) {
                    typedef typename ac_domain_vnodes::value_type v1;
                    remote_ac->insert(v1(p, vnode), true);
                } else {
                    typedef typename ac_domain_vnodes::value_type v1;
                    remote_ac->insert(v1(p, vnode), false);
                }
            }
            delete recv_mess;
        }

    // Calculate the widths of all the vnodes
    calcWidths();
}

//----------------------------------------------------------------------
// This differs from the previous ctor in that it allows preservation of
// global Vnode integer ID numbers associated with the input Vnodes. --tjw
template<unsigned Dim>
FieldLayout<Dim>::FieldLayout(const NDIndex<Dim>& domain,
			      const Vnode<Dim>* idx_begin,
			      const Vnode<Dim>* idx_end)
    : Domain(domain)
{

    // we have one more FieldLayout, indicate this
    //INCIPPLSTAT(incFieldLayouts);

    // for this kind of construction, we just take it that the user is
    // requesting all parallel axes
    for (unsigned int dl=0; dl < Dim; ++dl) RequestedLayout[dl] = PARALLEL;

    // Build Vnodes for each of the local domains.
    vnodesPerDirection_m = 0;
    int mynode = Ippl::Comm->myNode();
    for (const Vnode<Dim>* p=idx_begin; p!=idx_end; ++p) {
        typedef typename ac_id_vnodes::value_type v1;
        Local_ac.insert(v1(Unique::get(),
                           new Vnode<Dim>((*p).getDomain(),mynode,(*p).getVnode())));
    }

    // Everybody broadcasts their new local domains to everybody.
    // Build a message with the local domains and the id's for them.
    Message* bcast_mess = new Message;
    int count = idx_end - idx_begin;
    bcast_mess->put(count);
    // ::putMessage(*bcast_mess,idx_begin,idx_end);
    const Vnode<Dim>* indxmsg = idx_begin;
    while (indxmsg != idx_end)
        (indxmsg++)->putMessage(*bcast_mess);
    // Send it to everybody except yourself, and record the number sent.
    int tag = Ippl::Comm->next_tag( F_REPARTITION_BCAST_TAG , F_TAG_CYCLE );
    int node_count = Ippl::Comm->broadcast_others(bcast_mess,tag);
    // Create the container for the remote vnodes.
    ac_domain_vnodes *remote_ac = new ac_domain_vnodes( Domain );
    typedef typename ac_gc_domain_vnodes::value_type v1;
    Remotes_ac.insert(v1(gc0(),remote_ac) );

    // Loop until we receive a message from each other node.
    while ((--node_count)>=0)
        {
            // Receive a broadcast message from any node.
            int other_node = COMM_ANY_NODE;
            Message *recv_mess = Ippl::Comm->receive_block(other_node,tag);
            PAssert(recv_mess);
            // Extract the number of vnodes coming in.
            int count;
            recv_mess->get(count);
            // Now get the domains for the vnodes from the message.
            Vnode<Dim> p;
            while (count-- > 0) {
                p.getMessage(*recv_mess);
                Vnode<Dim> *vnode =
                    new Vnode<Dim>(p.getDomain(), other_node, p.getVnode());
                // For domains of extent 1 in all directions, must call
                // DomainMap::insert() with the noSplit flag set to true, to avoid an
                // assertion failure in Index::split. Also do this for domains of
                // extent ZERO. These are not supposed to happen in IPPL, but when the
                // do when something like BinaryRepartition creates them, accomodating
                // them here will prevent subsequent code hangs or other errors and
                // allow recovery such as rejection of the result of BinaryRepartition
                // (discard the FieldLayout with zero-size domains and go on with the
                // original layout unchanged). (tjw)
                //	if (p.getDomain().size() == 1) {
                if (p.getDomain().size() <= 1) {
                    typedef typename ac_domain_vnodes::value_type v1;
                    remote_ac->insert(v1(p.getDomain(), vnode), true);
                } else {
                    typedef typename ac_domain_vnodes::value_type v1;
                    remote_ac->insert(v1(p.getDomain(), vnode), false);
                }
            }
        }

    // Calculate the widths of all the vnodes
    calcWidths();
}


//----------------------------------------------------------------------
//
// Completely repartition this FieldLayout and all of the Fields
// defined on it.
//
template<unsigned Dim>
void
FieldLayout<Dim>::Repartition(const NDIndex<Dim>* idxBegin,
			      const NDIndex<Dim>* idxEnd)
{



    // Build a temporary FieldLayout to declare the temporary arrays.
    FieldLayout<Dim> tempLayout(Domain,idxBegin,idxEnd);

    // (TJW) Copy the vnodesPerDirection information:
    if (vnodesPerDirection_m != 0) {
        tempLayout.vnodesPerDirection_m = new unsigned[Dim];
        for (unsigned int d=0; d<Dim; d++) {
            tempLayout.vnodesPerDirection_m[d] = vnodesPerDirection_m[d];
        }
    }

    // Repartition each field.
    iterator_if p, endp=end_if();
    for (p=begin_if(); p!=endp; ++p) {
        FieldLayoutUser *user = (FieldLayoutUser *)((*p).second);
        user->Repartition(&tempLayout);
    }

    // Copy back the layout information.
    Local_ac = tempLayout.Local_ac;
    Remotes_ac = tempLayout.Remotes_ac;

    // Calculate the widths of all the vnodes
    calcWidths();

    //INCIPPLSTAT(incRepartitions);
}


//----------------------------------------------------------------------
//
// Completely repartition this FieldLayout and all of the Fields
// defined on it.
// This differs from the previous ctor in that it allows preservation of
// global Vnode integer ID numbers associated with the input Vnodes. --tjw
//
template<unsigned Dim>
void
FieldLayout<Dim>::Repartition(const Vnode<Dim>* idxBegin,
			      const Vnode<Dim>* idxEnd)
{



    // Build a temporary FieldLayout to declare the temporary arrays.
    FieldLayout<Dim> tempLayout(Domain,idxBegin,idxEnd);

    // (TJW) Copy the vnodesPerDirection information:
    if (vnodesPerDirection_m != 0) {
        tempLayout.vnodesPerDirection_m = new unsigned[Dim];
        for (unsigned int d=0; d<Dim; d++) {
            tempLayout.vnodesPerDirection_m[d] = vnodesPerDirection_m[d];
        }
    }

    // Repartition each field.
    iterator_if p, endp=end_if();
    for (p=begin_if(); p!=endp; ++p) {
        FieldLayoutUser *user = (FieldLayoutUser *)((*p).second);
        user->Repartition(&tempLayout);
    }

    // Copy back the layout information.
    Local_ac = tempLayout.Local_ac;
    Remotes_ac = tempLayout.Remotes_ac;

    // Calculate the widths of all the vnodes
    calcWidths();
}

//----------------------------------------------------------------------
//
// calculate the minimum vnode sizes in each dimension
template<unsigned Dim>
void FieldLayout<Dim>::calcWidths() {
    unsigned int d;

    // initialize widths first
    for (d=0; d < Dim; ++d)
        MinWidth[d] = getDomain()[d].length();

    // look for minimum width in local vnodes
    for  (const_iterator_iv v_i = begin_iv() ; v_i != end_iv(); ++v_i) {
        const NDIndex<Dim> &dom = (*v_i).second->getDomain();
        for (d=0; d < Dim; ++d) {
            if ((unsigned int) dom[d].length() < MinWidth[d])
                MinWidth[d] = dom[d].length();
        }
    }

    // look for minimum width in remove vnodes
    ac_domain_vnodes *v_ac = Remotes_ac[ gc0() ].get();
    for (iterator_dv dv_i = v_ac->begin(); dv_i != v_ac->end(); ++ dv_i) {
        const NDIndex<Dim> &dom = (*dv_i).first;
        for (d=0; d < Dim; ++d) {
            if ((unsigned int) dom[d].length() < MinWidth[d])
                MinWidth[d] = dom[d].length();
        }
    }

    //Inform dbgmsg("calcWidths", INFORM_ALL_NODES);
    //dbgmsg << "Calculated minimum widths for layout " << *this << endl;
    //for (d=0; d < Dim; ++d) {
    //  dbgmsg << " ==> MinWidth[" << d << "] = " << MinWidth[d];
    //  dbgmsg << " ... " << (getDistribution(d)==SERIAL?"serial":"parallel");
    //  dbgmsg << endl;
    //}
}


template<unsigned Dim>
NDIndex<Dim> FieldLayout<Dim>::getLocalNDIndex()
{
    NDIndex<Dim> theId;
    for (iterator_iv localv = begin_iv(); localv != end_iv(); ++localv) {
        Vnode<Dim> *vn = (*localv).second.get();
        if(vn->getNode() == Ippl::myNode())
            theId = vn->getDomain();
    }
    return theId;
}

//---------------------------------------------------------------------
// output
//---------------------------------------------------------------------
template<unsigned Dim>
void FieldLayout<Dim>::write(std::ostream& out) const
{


    int icount;

    // the whole domain, and the number of users
    out << "Domain = " << Domain << "\n";
    out << "FieldLayoutUsers = " << size_if() << "\n";

    // If applicable, vnodes per direction (tjw):
    if (vnodesPerDirection_m != 0) {
        out << "vnodesPerDirection_m[] =";
        for (unsigned int d=0; d<Dim; d++)
            out << " " << vnodesPerDirection_m[d];
        out << std::endl;
    }

    // iterate over the local vnodes and print them out.
    out << "Total number of vnodes = " << numVnodes() << std::endl;
    out << "Local Vnodes = " << Local_ac.size() << "\n";
    icount = 0;
    //tjw can now do operator<<() with whole Vnode objects
    //tjw  for(const_iterator_iv v_i = begin_iv() ; v_i != end_iv(); ++v_i)
    //tjw    out << " vnode " << icount++ <<" : "<< (*v_i).second->getDomain() << "\n";
    for(const_iterator_iv v_i = begin_iv() ; v_i != end_iv(); ++v_i)
        out << " vnode " << icount++ << ": " << *((*v_i).second) << std::endl;

    // iterate over the remote vnodes and print them out.
    ac_domain_vnodes *v_ac = Remotes_ac[ gc0() ].get();
    out << "Remote Vnodes = " << v_ac->size() << "\n";
    icount = 0;
    for (iterator_dv dv_i = v_ac->begin(); dv_i != v_ac->end(); ++ dv_i)
        out << " vnode " << icount++ << " : " << *((*dv_i).second) << std::endl;
}

// not in IPPL 2
//---------------------------------------------------------------------
template<unsigned Dim>
void FieldLayout<Dim>::new_gc_layout(const GuardCellSizes<Dim>& gc)
{
    // Build the guarded domain.
    NDIndex<Dim> guarded_domain( AddGuardCells(Domain,gc) );
    // Build a container for vnodes in that domain.
    ac_domain_vnodes *gr = new ac_domain_vnodes(guarded_domain);
    // Record pointer to that container using gc as the key.
    typedef typename ac_gc_domain_vnodes::value_type v1;
    Remotes_ac.insert(v1(gc,gr) );
    // Get the container of vnodes stored w/o guard cells.
    ac_domain_vnodes &v0 = *Remotes_ac[ gc0() ];
    // Loop over all the remote vnodes.
    for (iterator_dv v_i = v0.begin(); v_i != v0.end(); ++v_i)
        {
            // Build the domain for this vnode with gc guard cells.
            NDIndex<Dim> domain(AddGuardCells((*v_i).first,gc));
            // Record pointer to this vnode with these guard cells.
            typedef typename ac_domain_vnodes::value_type v2;
            gr->insert(v2(domain,(*v_i).second) );
        }
}

//----------------------------------------------------------------------
//
// Tell the FieldLayout that a FieldLayoutUser is using it
template<unsigned Dim>
void
FieldLayout<Dim>::checkin(FieldLayoutUser& f,
			  const GuardCellSizes<Dim>& gc)
{
    checkinUser(f);
    iterator_gdv guarded = Remotes_ac.find(gc);
    if ( guarded == Remotes_ac.end() )
        new_gc_layout(gc);
}

//----------------------------------------------------------------------
//
// Tell the FieldLayout that a Field is no longer using it.
template<unsigned Dim>
void
FieldLayout<Dim>::checkout(FieldLayoutUser& f)
{
    checkoutUser(f);
}

