//
// Class Format
// Format class to allow serializing message objects into plain buffers
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
#include "Message/Format.h"

Format::Format(Message *msg)
{
    items = msg->size();
    size = 0;

    format_array.resize(2*items);
    for (unsigned int i=0; i<items; ++i)
    {
        Message::MsgItem &msgitem = msg->item(i);
        format_array[2*i+0] = msgitem.numElems();
        format_array[2*i+1] = msgitem.numBytes();
        size += format_array[2*i+1];
    }
}

void Format::print()
{
	std::cout << "size: " << size << std::endl;
	for (unsigned int i=0; i<items; ++i)
    {
		std::cout << "entry " << i << ": " << format_array[2*i+0]
			<< " elements " << format_array[2*i+1] << " bytes\n";
    }
}