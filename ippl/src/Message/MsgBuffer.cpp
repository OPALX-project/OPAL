//
// Class MsgBuffer
// MsgBuffer class to allow serializing message objects into plain buffers
// to send directly with mpi calls or similar means
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
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
#include "Message/MsgBuffer.h"

#include "Message/Format.h"

MsgBuffer::MsgBuffer(Format *f, int count, int offset)
        : format(f), writepos(0), readpos(0)
{
    datasize = count*format->getSize();
    data.resize(datasize+offset);
}


MsgBuffer::MsgBuffer(Format *f, char *buf, int size)
        : format(f), writepos(0), readpos(0)
{
    datasize = size;
    data.resize(datasize);
    std::copy(buf, buf+size, data.begin());
    delete[] buf;
}

MsgBuffer::~MsgBuffer()
{
}

bool MsgBuffer::add(Message *msg)
{
    unsigned int items = msg->size();

    //check for full storage or message size mismatch
    if (writepos == datasize || items != format->getItemCount())
        return false;

    int pos = writepos;
    for (unsigned int i=0; i<items; ++i)
    {
        Message::MsgItem &msgitem = msg->item(i);

        //check for format mismatch
        if (format->getItemElems(i) != msgitem.numElems() ||
                format->getItemBytes(i) != msgitem.numBytes())
            return false;

        //actually copy to buffer
        std::memcpy(data.data()+pos, msgitem.data(), format->getItemBytes(i));
        pos += format->getItemBytes(i);
    }

    writepos = pos;
    return true;
}

Message* MsgBuffer::get()
{
    if (readpos > datasize - format->getSize())
        return 0;

    unsigned int items = format->getItemCount();
    Message *msg = new Message(items);

    //get all the items according to format and add them to the message
    for (unsigned int j = 0; j < items; j++)
    {
        unsigned int bytesize = format->getItemBytes(j);
        unsigned int elements = format->getItemElems(j);

        msg->setCopy(false);
        msg->setDelete(false);
        msg->putmsg(data.data()+readpos, bytesize/elements, elements);
        readpos += bytesize;
    }

    return msg;
}