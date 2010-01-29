/*
  NAME:
     calc_qstarij (NOT USED ANYMORE)
  PURPOSE:
     calculate the qstar_sum_factor for the partial E step
  CALLING SEQUENCE:
     calc_qstarij(double * qstarij, gsl_matrix * qij, int partial_indx[3])
  INPUT:
     qij          - qij vector of old posterior loglikelihoods (given as logs)
     partial_indx - indices of the partial EM gaussians
  OUTPUT:
     qstarij - sum_j=1,2,3 of q*ij
  REVISION HISTORY:
     2008-09-22 - Written Bovy
*/
#include <math.h>
#include <float.h>
#include <gsl/gsl_matrix.h>
#include <proj_gauss_mixtures.h>

void calc_qstarij(double * qstarij, gsl_matrix * qij, int partial_indx[3]){
  int ii;
  for (ii= 0; ii != qij->size1; ++ii)
    *(qstarij++) = logsum(qij,ii,true,true,partial_indx);
      /*DBL_MIN+exp(gsl_matrix_get(qij,ii,partial_indx[0])) +
      exp(gsl_matrix_get(qij,ii,partial_indx[1])) +
      exp(gsl_matrix_get(qij,ii,partial_indx[2]));
      */
  qstarij -= qij->size1;

  return;
}
