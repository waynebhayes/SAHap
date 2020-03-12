#include "utils.hpp"

/*
function LogPoissonPMF(l,k, pr,r,i){
  pr=0;r=-l;
  for(i=k;i>0;i--){
    pr=r;r+=log(l/i)
  }
  return r
}
*/
double log_poisson_pmf(double l, unsigned k) {
  double r = -l;

  for (auto i = k; i > 0; i--) {
    r += log(l/i);
  }

  return r;
}

double log_poisson_1_cdf_l = 0;
double log_poisson_1_cdf_memo[1000];

/*
function LogPoisson1_CDF(l,k, i,sum,psum){
  pmax=2;max=-1e30;
  for(i=k;pmax!=max;i++){
    pmax=max;max=MAX(max,LogPoissonPMF(l,i))
  };
  if(max==1 && k<l) return 0; # this means the numbers are so big the sum got zero but we got less than expected.
  else return max/.894
}
*/
double log_poisson_1_cdf(double l, unsigned k) {
  /*
  if (log_poisson_1_cdf_l != l) {
    fill(log_poisson_1_cdf_memo, log_poisson_1_cdf_memo + 1000, 10);
    log_poisson_1_cdf_l = l;
  }

  if (k < 1000 && log_poisson_1_cdf_memo[k] < 0) {
    return log_poisson_1_cdf_memo[k];
  }
  */
  // cout << "log_poisson_1_cdf(" << l << ", " << k << ")" << endl;

  double pmax = 2, max = -1e30;
  for (unsigned i = k; pmax != max; ++i) {
    pmax = max;
    double logpmf = log_poisson_pmf(l, i);
    if (logpmf > max) {
      max = logpmf;
    }
  }

  double r = (max == 1 && k < l) ? 0 : max / .894;

  /*
  if (k < 1000) {
    log_poisson_1_cdf_memo[k] = r;
  }

  cout << "k = " << k << ", r = " << r << endl;
  */

  return r;
}