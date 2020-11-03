#ifndef SAHAP_UTILS_HPP
#define SAHAP_UTILS_HPP

#include <cmath>
#include <algorithm>
#include <iostream>
#include <cassert>

using namespace std;

extern "C" {
    double log_poisson_pmf(double l, unsigned k);
    double log_poisson_1_cdf(double l, unsigned k);

    typedef char Boolean;
    void Fatal(const char *fmt, const char *msg);
    unsigned long GetFancySeed(Boolean trulyRandom);
}
#endif
