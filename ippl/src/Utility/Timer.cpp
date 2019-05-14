#include "Timer.h"

Timer::Timer() {
    this->clear();
}


void Timer::clear() {
    wall_m = user_m = sys_m = 0.0;
}


void Timer::start() {
    timer_m.start();
}


void Timer::stop() {
    timer_m.stop();
    
    boost::timer::cpu_times elapsed = timer_m.elapsed();
    
    wall_m += elapsed.wall;
    user_m += elapsed.user;
    sys_m  += elapsed.system;
}


double Timer::clock_time() {
    return wall_m * 1.0e-9;
}


double Timer::user_time() {
    return user_m * 1.0e-9;
}


double Timer::system_time() {
    return sys_m * 1.0e-9;
}


double Timer::cpu_time() {
    return (user_m + sys_m) * 1.0e-9;
}
