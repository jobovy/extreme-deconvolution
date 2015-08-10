#include <proj_gauss_mixtures_IDL.h>
#include <stdio.h>

char * bool2char(int * a, int K)
{
	char x[K];

	for (int kk = 0; kk != K; ++kk) {
		x[kk] = (char)*(a + kk);
	}
	char * px = x;
	return px;
}


void proj_gauss_mixtures_R(double * ydata, double * ycovar,
                           double * projection, double * logweights,
                           int * N, int * dy,
                           double * amp, double * xmean,
                           double * xcovar, int * d, int * K,
                           int * fixamp, int * fixmean,
                           int * fixcovar,
                           double * avgloglikedata, double * tol,
                           int * maxiter, int * likeonly, double * w,
                           char * logfilename, int * slen, int * splitnmerge,
                           char * convlogfilename, int * convloglen,
                           int * noprojection, int * diagerrors,
                           int * noweights)
{
	char * char_fixamp = bool2char(fixamp, *K);
	char * char_fixmean = bool2char(fixmean, *K);
	char * char_fixcovar = bool2char(fixcovar, *K);

	proj_gauss_mixtures_IDL(ydata, ycovar, projection, logweights,
		*N, *dy, amp, xmean, xcovar, *d, *K, char_fixamp, char_fixmean,
		char_fixcovar, avgloglikedata, *tol, *maxiter, (char)(*likeonly), *w,
		logfilename, *slen, *splitnmerge, convlogfilename, *convloglen,
		(char)(*noprojection), (char)(*diagerrors), (char)(*noweights));
	/* int res_size = sizeof(avgloglikedata)/sizeof(double); */
	/* for (int kk = 0; kk != res_size; ++kk) printf("%f\n", avgloglikedata[kk]); */
}


