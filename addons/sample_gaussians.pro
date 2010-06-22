;+
;   NAME:
;      sample_gaussians
;
;   PURPOSE:
;      sample a distribution that is the some of gaussians
;
;   CALLING SEQUENCE:
;      sample_gaussians, nsamples=nsamples, mean=mean, covar=covar,
;      amp=amp, sample=sample
;
;   INPUT:
;      nsamples - number of points desired
;      mean     - array of means [d,M]
;      covar    - array of covariances [d,d,M]
;      amp      - array of amplitudes [M]
;
;   OPTIONAL INPUTS:
;      seed    - seed to use to initialize the random series
;   OUTPUT:
;      sample - array of sampled points [d,N]
;
;   REVISION HISTORY:
;      2008-07-30 - Written Bovy 
;-
PRO SAMPLE_GAUSSIANS, nsamples=nsamples, mean=mean, covar=covar, amp=amp, $
                      sample=sample, seed=seed

ON_ERROR, 2

xmean=mean
xcovar=covar
xamp= amp

;;dimension+ngauss
d= n_elements(xmean)/n_elements(xamp)
ngauss= n_elements(xamp)

IF ~keyword_set(nsamples) THEN nsamples= 1
IF ~keyword_set(seed) THEN seed= 1L

IF d EQ 1 THEN BEGIN
    ;;Accumulate amp
    xamp = TOTAL(xamp,/cumulative,/double)
    ;;Then loop over nsamples
    sample= dblarr(nsamples)
    FOR ii= 0L, nsamples-1 DO BEGIN
        gauss= RANDOMU(seed,/double)
        jj= 0L
        WHILE (gauss GT xamp[jj]) DO jj++
        sample[ii]= sqrt(xcovar[jj])*RANDOMN(seed,/double) + xmean[jj]
    ENDFOR
    RETURN
ENDIF


;;First diagonalize all covariances
var= dblarr(d,ngauss)
FOR gg=0L, ngauss-1 DO BEGIN
    ;;IDL sux and therefore we must symmetrize xcovar
    xcovar[*,*,gg]= .5D*(xcovar[*,*,gg]+transpose(xcovar[*,*,gg]))
    var[*,gg]= eigenql(xcovar[*,*,gg],/double,eigenvectors=ev)
    xcovar[*,*,gg]= ev
ENDFOR

;;Accumulate amp
xamp = TOTAL(xamp,/cumulative,/double)

;;Then loop over nsamples
sample= dblarr(d,nsamples)
FOR ii= 0L, nsamples-1 DO BEGIN
    gauss= RANDOMU(seed,/double)
    jj= 0L
    WHILE (gauss GT xamp[jj]) DO jj++
    tmpsample= dblarr(d)
    FOR kk=0L, d-1 DO tmpsample[kk]= sqrt(var[kk,jj])*RANDOMN(seed,/double)
    tmpsample= transpose(xcovar[*,*,jj])##tmpsample
    sample[*,ii]= tmpsample+xmean[*,jj]
ENDFOR

END
