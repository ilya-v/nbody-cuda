#ifndef STATUS_HPP_
#define STATUS_HPP_

#include <string>

struct status_t {
    unsigned step;
    double t, dt, ep, ek, etot, I;
    double machine_time;
    long rng;

    std::string header() const;
    std::string line() const;
};

#endif // !STATUS_HPP_

