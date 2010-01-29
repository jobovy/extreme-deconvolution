/*
  NAME:
     normalize_row
  PURPOSE:
     normalize a row (or column) of a matrix given as logs
  CALLING SEQUENCE:
     normalize_row(gsl_matrix * q, int row,bool partial, int partial_indx[3])
  INPUT:
     q            - matrix
     row          - row to be normalized
     isrow        - is it a row or a column
  OUTPUT:
     normalization factor (i.e. logsum)
  REVISION HISTORY:
     2008-09-21 - Written Bovy
*/
#include <gsl/gsl_matrix.h>
#include <proj_gauss_mixtures.h>

double normalize_row(gsl_matrix * q, int row, bool isrow){
  double loglike;
  if (isrow)
    loglike = logsum(q,row,true);
  else
    loglike = logsum(q,row,false);

  int dd;
  if (isrow)
    for (dd = 0; dd != q->size2; ++dd)
      gsl_matrix_set(q,row,dd,gsl_matrix_get(q,row,dd)-loglike);
  else
    for (dd = 0; dd != q->size1; ++dd)
      gsl_matrix_set(q,dd,row,gsl_matrix_get(q,dd,row)-loglike);
    

  return loglike;
}
