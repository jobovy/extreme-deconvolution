;+
; NAME:
;    exd_jackerr_example
; PURPOSE:
;    example of the use of exd_jackerr
; CALLING SEQUENCE:
; INPUT:
; OUTPUT:
; HISTORY:
;    2013-01-05 - Written - Bovy (IAS)
;-
PRO EXD_JACKERR_EXAMPLE

;;sample some data from 3 Gaussians
real_amp= dblarr(3)
real_amp[0]= 0.3
real_amp[1]= 0.4
real_amp[2]= 0.3
real_mean= dblarr(1,3)
real_mean[0,0]= -2.
real_mean[0,1]= 2.
real_mean[0,2]= 0.
real_covar= dblarr(1,1,3)
real_covar[0,0,0]= 1.
real_covar[0,0,1]= 1.
real_covar[0,0,2]= 1.
nsamples= 1000L
seed= 1L
sample_gaussians, nsamples=nsamples, mean=real_mean, covar=real_covar, $
  amp=real_amp, sample=sample, seed=seed
sample= REFORM(sample,1,nsamples)
sample_covar= dblarr(1,nsamples)

for kk=1L, 5 do begin
    ;;Fit with 1 Gaussian using exd_jackerr
    ngauss= kk
    njack= 30 ;nsamples
    fit_amp= dblarr(ngauss)
    fit_mean= dblarr(1,ngauss)
    fit_covar= dblarr(1,1,ngauss)
    FOR ii=0L, ngauss-1 DO BEGIN
        fit_amp[ii]= 1./double(ngauss)
        fit_mean[0,ii]= randomn(seed,1)
        fit_covar[0,0,ii]= 4.
    ENDFOR
    projected_gauss_mixtures_c, ngauss, sample, sample_covar, $
      fit_amp, fit_mean, fit_covar, /quiet
    exd_jackerr, ngauss, sample, sample_covar, $
      fit_amp, fit_mean, fit_covar, $
      fit_amp_err, fit_mean_err, fit_covar_err, $
      jacklogl=jacklogl, njack=njack
    print, "jack logl: ", jacklogl
    print, "amplitude and  error: ", fit_amp, fit_amp_err
endfor

END
