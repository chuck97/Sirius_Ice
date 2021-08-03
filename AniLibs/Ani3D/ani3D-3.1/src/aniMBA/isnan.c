//to call floating-point check from the Fortran routine
#include <math.h>

void fpcheck_( double *a, int *flag ) 
{ 
   *flag = isnan(*a) && !isinf(*a); 
}

