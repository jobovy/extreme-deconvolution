/*
  NAME:
     bovy_isfin
  PURPOSE:
     Is the given value finite? (i.e. -DBL_MAX <= x <= DBL_MAX
  CALLING SEQUENCE:
     finite = bovy_isfin(double x)
  INPUT:
     x  -  double
  OUTPUT:
     true if x is finite
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <stdbool.h>
#include <float.h>

inline bool bovy_isfin(double x){
  return (x > DBL_MAX || x < -DBL_MAX) ? false : true;
}
