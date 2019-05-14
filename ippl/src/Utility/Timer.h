#ifndef TIMER_H
#define TIMER_H

#include <boost/timer/timer.hpp>

// https://www.boost.org/doc/libs/1_70_0/libs/timer/doc/cpu_timers.html
class Timer
{
public:
    
    Timer();
    
    void clear();               // Set all accumulated times to 0
    void start();               // Start timer
    void stop();                // Stop timer
    
    double clock_time();        // Report clock time accumulated in seconds
    double user_time();         // Report user time accumlated in seconds
    double system_time();       // Report system time accumulated in seconds
    double cpu_time();          // Report total cpu_time which is just user_time + system_time
    
private:
    double wall_m;
    double user_m;
    double sys_m;
    boost::timer::cpu_timer timer_m;
};

#endif
