#include <unistd.h>             /* sysconf(_SC_CLK_TCK) */
#include <sys/times.h>          /* struct tms, times() */


#define SECS_PER_CLOCK (1./sysconf(_SC_CLK_TCK))


//void seconds_( double *usecs, double *ssecs ) 
void seconds_( double *usecs ) 
{
   struct tms t;
   clock_t utime;
   clock_t stime;

   times(&t);
   utime = t.tms_utime;
   stime = t.tms_stime;

   *usecs = utime * SECS_PER_CLOCK;
   // *ssecs = stime * SECS_PER_CLOCK;
}

