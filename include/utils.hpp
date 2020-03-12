#ifndef SAHAP_UTILS_HPP
#define SAHAP_UTILS_HPP

#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

double log_poisson_pmf(double l, unsigned k);
double log_poisson_1_cdf(double l, unsigned k);

#endif