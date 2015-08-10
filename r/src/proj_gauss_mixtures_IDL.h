#ifndef __PROJ_GAUSS_MIXTURES_IDL__
#define __PROJ_GAUSS_MIXTURES_IDL__
int proj_gauss_mixtures_IDL(double * ydata, double * ycovar, 
			    double * projection, double * logweights,
			    int N, int dy, 
			    double * amp, double * xmean, 
			    double * xcovar, int d, int K, 
			    char * fixamp, char * fixmean, 
			    char * fixcovar, 
			    double * avgloglikedata, double tol, 
			    int maxiter, char likeonly, double w, 
			    char * logfilename, int slen, int splitnmerge,
			    char * convlogfilename, int convloglen,
			    char noprojection,char diagerrors,
			    char noweights);
#endif /* proj_gauss_mixtures_IDL.h */
