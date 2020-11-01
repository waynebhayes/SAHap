#include "utils.hpp"
#include <unistd.h>
#include <string.h>

#if defined(__WATCOMC__) || defined(__MINGW32__) || defined(__CYGWIN__)
#define POPEN 0
#else
#define POPEN 1
#endif

extern "C" {

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

void Fatal(const char *fmt, const char *msg) { fprintf(stderr, fmt, msg); exit(1);}

FILE *Popen(const char *cmd, const char *mode) {
#if POPEN
    return popen(cmd, mode);
#else
    char CMD[strlen(cmd) + 100], TMP[100];
    sprintf(TMP, "/tmp/Popen.%d",getpid());
    sprintf(CMD, "(%s) > %s; (sleep 10; /bin/rm %s) &",cmd, TMP,TMP);
    system(CMD);
    return fopen(TMP,mode);
#endif
}
void  Pclose(FILE *fp) {
#if POPEN
    pclose(fp);
#else
    fclose(fp);
#endif
}

/* From Libwayne:
** Try to compute a seed that will be different for all processes even if they're all started at
** the same time, on the same or different servers. We use the host's IPv4 address, the time
** (to the nearest second), the process ID, and the parent process ID. The only gotcha is that
** if you call this twice within the same second within the same process, the result will be the
** same. But since you should *never* seed twice within the same code, that's your problem.
** (This problem can be offset by setting "trulyRandom" to true.)
*/
unsigned long GetFancySeed(Boolean trulyRandom)
{
    unsigned long seed = 0;
    const char *cmd = "hostname -i | awk '{for(i=1;i<=NF;i++)if(match($i,\"^[0-9]*\\\\.[0-9]*\\\\.[0-9]*\\\\.[0-9]*$\")){IP=$i;exit}}END{if(!IP)IP=\"127.0.0.1\"; print IP}'";

    FILE *fp=Popen(cmd,"r");
    int i, ip[4], host_ip=0;
    if(4!=fscanf(fp," %d.%d.%d.%d ", ip, ip+1, ip+2, ip+3)) Fatal("Attempt to get IPv4 address failed:\n%s\n",cmd);
    Pclose(fp);
    for(i=0;i<4;i++) host_ip = 256*host_ip + ip[i];
    unsigned long dev_random=0;
    if(trulyRandom) {
	fp = fopen("/dev/urandom","r");
	if(!fp) fp = fopen("/dev/random","r");
	if(fp){
	    assert(1 == fread(&dev_random, sizeof(dev_random),1, fp));
	    fclose(fp);
	}
	else
	{
	    // Use a bunch of hard-to-predict commands with nondeterministic output, then md5sum it.
	    fp = Popen("(who; w; uptime; ipconfig /all; ifconfig -a) 2>/dev/null | md5sum","r");
	    unsigned char c;
	    while((c=fgetc(fp)) > 0 && c != ' ' && c != '\t') {
		unsigned int hex;
		if(isdigit(c)) hex = (c-'0');
		else {
		    c = tolower(c);
		    assert(islower(c) && c >= 'a' && c <= 'f');
		    hex = 10+(c-('a'-1));
		}
		dev_random = 16*dev_random + hex;
	    }
	    Pclose(fp);
	}
    }
    seed = host_ip + time(0) + getppid() + getpid() + dev_random;
#if 0
    fprintf(stderr,"%s\n",cmd);
    fprintf(stderr,"%d.%d.%d.%d\n",ip[0],ip[1],ip[2],ip[3]);
    fprintf(stderr,"seed is %ud\n",seed);
#endif
    return seed;
}

} // extern "C"
