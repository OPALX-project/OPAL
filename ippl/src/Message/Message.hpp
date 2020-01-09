// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Message/Message.h"
#include "Utility/Inform.h"
#include "Utility/IpplInfo.h"


#include <iterator>
// Rogue Wave STL for PGI compiler doesn't have iterator_traits!!

// By now they have 2002 andreas
#ifdef OLD_IPPL_PGI
namespace std
{

template <class _Iterator>
struct iterator_traits
{
    typedef typename _Iterator::iterator_category iterator_category;
    typedef typename _Iterator::value_type        value_type;
    typedef typename _Iterator::difference_type   difference_type;
    typedef typename _Iterator::pointer           pointer;
    typedef typename _Iterator::reference         reference;
};

template <class _Tp>
struct iterator_traits<_Tp*>
{
    typedef random_access_iterator_tag iterator_category;
    typedef _Tp                         value_type;
    typedef ptrdiff_t                   difference_type;
    typedef _Tp*                        pointer;
    typedef _Tp&                        reference;
};

template <class _Tp>
struct iterator_traits<const _Tp*>
{
    typedef random_access_iterator_tag iterator_category;
    typedef _Tp                         value_type;
    typedef ptrdiff_t                   difference_type;
    typedef const _Tp*                  pointer;
    typedef const _Tp&                  reference;
};

}
#endif // IPPL_PGI



#include <new>
#include <memory>
#include <cstdlib>
#include <cstddef>

// MWERKS: moved macro definitions for MessageTypeIntrinsic to Message.h

// MWERKS: moved PutSingleItem definitions to Message.h

// definitions for non-templated static member functions of PutSingleItem
// can remain here.  JCC

// specialization to a non-built-in type, which is never assumed to be a ptr

template <class T>
Message&
PutSingleItem<T, false, false>::put(Message& msg, const T& cvalue)
{
    T& value = const_cast<T&>(cvalue);
    value.putMessage(msg);
    return msg;
}

template <class T>
Message&
PutSingleItem<T, false, false>::get(Message& msg, T& value)
{
    value.getMessage(msg);
    return msg;
}

template <class T>
Message&
PutSingleItem<T, false, false>::put(Message& m, T beg, T end)
{
    // get the value type using iterator traits
    typedef typename std::iterator_traits<T>::value_type T2;

    // make sure this is not empty ... if so, just create an empty item
    if (beg == end)
    {
        m.putmsg(0, sizeof(T2), 0);
    }
    else
    {
        // find the number of elements in this iterator range
        unsigned int d = 0;
        T f;
        for (f = beg; f != end; ++d, ++f);

        // if DoCopy is false, we must assume the iterators are pointers
        // to simple data types, and so we can just cast the pointer to
        // void and store it
        if (!m.willCopy())
        {
            m.putmsg((void*) (&*beg), sizeof(T2), d);
        }
        else
        {
            // make a copy ourselves
            // use malloc to get block, placement new for each element
            T2* cpydata = static_cast<T2*>( malloc(sizeof(T2) * d) );
            T2* cpy = cpydata;
            T i;
            for (i = beg; i != end; ++i, ++cpy)
                new (cpy) T2(*i);

            // put data into this message
            m.setCopy(false);
            m.setDelete(true);
            m.putmsg( (void*) cpydata, sizeof(T2), d );
        }
    }
    return m;
}

template <class T>
Message&
PutSingleItem<T, false, false>::put(Message& m,
                                    const std::vector<size_t>& indices, T beg)
{
    // get the value type using iterator traits
    typedef typename std::iterator_traits<T>::value_type T2;

    // make sure this is not empty ... if so, just create an empty item
    if (indices.begin() == indices.end())
    {
        m.putmsg(0, sizeof(T2), 0);
    }
    else
    {
        // find the number of elements in this iterator range
        std::vector<size_t>::size_type d = indices.size();

        // make a copy of the data ourselves
        // use malloc to get block, placement new for each element
        T2* cpydata = static_cast<T2*>( malloc(sizeof(T2) * d) );
        T2* cpy = cpydata;
        std::vector<size_t>::const_iterator i, iend = indices.end();
        for (i = indices.begin(); i != iend; ++i, ++cpy)
            new (cpy) T2(beg[*i]);

        // put data into this message
        m.setCopy(false);
        m.setDelete(true);
        m.putmsg( (void*) cpydata, sizeof(T2), d );
    }
    return m;
}

template <class T>
Message&
PutSingleItem<T, false, false>::get_iter(Message& m, T o)
{
    // get the value type using iterator traits
    typedef typename std::iterator_traits<T>::value_type T2;

    // check to see if there is an item
    if ( m.empty() )
    {
        ERRORMSG("get_iter(): no more items in Message" << endl);
    }
    else
    {
        // get the next MsgItem off the top of the list
        Message::MsgItem& mitem = m.item(0);

        // copy the data pointer
        T2* data = static_cast<T2*>( mitem.data() );
        int i;
        for (i = mitem.numElems(); i > 0; i--)
            *o++ = *data++;

        // delete this MsgItem
        m.get();
    }
    return m;
}


// specialization to a built-in type that is not a pointer

template <class T>
Message&
PutSingleItem<T, true, false>::put(Message& msg, const T& value)
{
    T& ncval = const_cast<T&>(value);
    return msg.putmsg( (void*) &ncval, sizeof(T) );
}

template <class T>
Message&
PutSingleItem<T, true, false>::get(Message& msg, T& value)
{
    return msg.getmsg( (void*) &value );
}


// specialization to a pointer to a built-in type. In this class, we
// know that 'T' is a pointer type.


template <class T>
Message&
PutSingleItem<T, true, true>::put(Message& msg, T beg, T end)
{
    typedef typename std::iterator_traits<T>::value_type value_type;
    return msg.putmsg( (void*) beg,
                       sizeof(value_type),
                       (end - beg) );
}

template <class T>
Message&
PutSingleItem<T, true, true>::get(Message& msg, T ptr)
{
    return msg.getmsg( (void*) ptr );
}

template <class T>
Message&
PutSingleItem<T, true, true>::put(Message& m,
                                  const std::vector<size_t>& indices, T beg)
{
    return PutSingleItem<T, false, false>::put(m, indices, beg);
}

template <class T>
Message&
PutSingleItem<T, true, true>::get_iter(Message& m, T o)
{
    return PutSingleItem<T, false, false>::get_iter(m, o);
}
