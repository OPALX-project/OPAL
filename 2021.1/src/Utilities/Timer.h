#ifndef OPAL_Timer_HH
#define OPAL_Timer_HH 1

// ------------------------------------------------------------------------
// $RCSfile: Timer.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Timer
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:48 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include <ctime>
#include <string>


// Class Timer
// ------------------------------------------------------------------------
/// Timer class.
//  Encapsulates the UNIX wall time functions.
//  Get Calendar date or time.

namespace OPALTimer {

    class Timer {

    public:

        /// Constructor.
        //  Store the system clock time.
        Timer();

        ~Timer();

        /// Return date.
        std::string date() const;

        /// Return time.
        std::string time() const;

    private:

        // Not implemented.
        Timer(const Timer &);
        void operator=(const Timer &);

        time_t timer;
    };

}
#endif // OPAL_Timer_HH
