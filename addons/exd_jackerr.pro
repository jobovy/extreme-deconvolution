;+
; NAME:
;   exd_jackerr
; PURPOSE:
;   use jackknife to estimate uncertainties on the best-fit parameters
;   of the Gaussian mixture
; CALLING SEQUENCE:
;   exd_jackerr, ngauss, ydata, ycovar, $
;                amp, xmean, xcovar, $
;                amperr, xmeanerr, xcovarerr, $
;                projection=projection, $
;                fixamp=fixamp,fixmean=fixmean, $
;                fixcovar=fixcovar, $
;                tol=tol, maxiter=maxiter, $
;                w=w, $
;                quiet=quiet, splitnmerge=splitnmerge, maxsnm=maxsnm, $
;                weight=weight, logweight=logweight
; INPUTS:
;   ngauss     - number of desired gaussians to fit
;   ydata      - [ndimy,ndata] observed velocities
;   ycovar     - [ndimy(,ndimy),ndata] observed velocities' errors-squared 
;                (if [ndimy,ndata] then the errors are assumed to be 
;                uncorrelated)
;   amp        - [ngauss] list of relative amplitudes (best-fit)
;   xmean      - [ndimx,ngauss] list of initial guesses for the means
;                (best-fit)
;   xcovar     - [ndimx,ndimx,ngauss] list of initial guesses for the
;                 covar (best-fit)
;
; OPTIONAL INPUTS:
;   njack      - number of jackknife subsamples to use (default: ndata)
;   projection - [ndimx, ndimy, ndata] non-square matrices
;                implementing a projection from the model space to the
;                data space for each data point
;   weight     - weights to be applied to the data points [ndimy]
;   fixamp     - [ngauss] list of integers: 0 to update amp, 1 not to
;   fixmean    - [ngauss] list of integers: 0 to update xmean, 1 not to
;   fixcovar   - [ngauss] list of integers: 0 to update xcovar, 1 not
;                to
;   tol        - tolerance (convergence iff difference in avgloglike <
;                tol)
;   maxiter    - maximum number of iterations
;   w          - constrain any variance to be larger than this
;                value
;   splitnmerge- depth to go down the splitnmerge path
;
; KEYWORDS:
;   maxsnm     - use the maximum number of split 'n' merge steps,
;                K*(K-1)*(K-2)/2
;   logweight  - indicates that the weights in weight are actually ln(weight)
;
; OUTPUTS:
;   xmeanerr, xcovarerr and amperr
; REVISION HISTORY:
;   2010-06-11  - started on one-d case- Bovy (NYU)
;-
PRO EXD_JACKERR, ngauss, ydata, ycovar, $
                 amp, xmean, xcovar, $
                 amperr, xmeanerr, xcovarerr, $
                 njack=njack, _EXTRA=_EXDKEYS

ndimy=(size(ydata,/dimensions))[0]
ndata=n_elements(ydata)/ndimy
ndimx=n_elements(xmean)/n_elements(amp)
IF n_elements(ycovar) EQ n_elements(ydata) THEN BEGIN
    diagerrors= 1B
    ycovar= reform(ycovar,ndimx,ndata)
ENDIF ELSE BEGIN
    diagerrors= 0B
    ycovar= reform(ycovar,ndimx,ndimx,ndata)
ENDELSE

IF ~keyword_set(njack) THEN njack= ndata
IF njack NE ndata THEN BEGIN
    ndatasub= floor(double(ndata)/njack)
    njack= ceil(double(ndata)/ndatasub)
ENDIF

;;Run ExD for each jackknife subsample and collect the results
initamp= amp
initxmean= xmean
initxcovar= xcovar

amps= dblarr(ngauss,njack)
xmeans= dblarr(ndimx,ngauss,njack)
xcovars= dblarr(ndimx,ndimx,ngauss,njack)
FOR ii=0L, njack-1 DO BEGIN
    IF njack EQ ndata THEN BEGIN
        ;;ydata
        IF ii EQ 0 THEN BEGIN
            thisydata= ydata[*,1:ndata-1]
        ENDIF ELSE IF ii EQ ndata-1 THEN BEGIN
            thisydata= ydata[*,0:ndata-2]
        ENDIF ELSE BEGIN
            thisydata= transpose([transpose(ydata[*,0:ii-1]),$
                                  transpose(ydata[*,ii+1:ndata-1])])
        ENDELSE
        thisydata= reform(thisydata,ndimy,ndata-1)
        ;;ycovar
        IF diagerrors THEN BEGIN
            IF ii EQ 0 THEN BEGIN
                thisycovar= ycovar[*,1:ndata-1]
            ENDIF ELSE IF ii EQ ndata-1 THEN BEGIN
                thisycovar= ycovar[*,0:ndata-2]
            ENDIF ELSE BEGIN
                thisycovar= transpose([transpose(ycovar[*,0:ii-1]),$
                                       transpose(ycovar[*,ii+1:ndata-1])])
            ENDELSE
            thisycovar= reform(thisycovar,ndimy,ndata-1)
        ENDIF ELSE BEGIN 
            IF ii EQ 0 THEN BEGIN
                thisycovar= ycovar[*,*,1:ndata-1]
            ENDIF ELSE IF ii EQ ndata-1 THEN BEGIN
                thisycovar= ycovar[*,*,0:ndata-2]
            ENDIF ELSE BEGIN
                thisycovar= transpose([transpose(ycovar[*,*,0:ii-1]),$
                                       transpose(ycovar[*,*,ii+1:ndata-1])])
            ENDELSE
            thisycovar= reform(thisycovar,ndimy,ndimy,ndata-1)
        ENDELSE
    ENDIF ELSE BEGIN
        print, "njack NE ndata not currently implemented"
    ENDELSE
    thisamp= initamp
    thisxmean= initxmean
    thisxcovar= initxcovar
    ;;Run ExD on the jackknife sample
    print, format = '("Working on jackknife sample ",i7," of ",i7,a1,$)', $
      ii+1,njack,string(13B)
    avgloglikedata= 0D0
    projected_gauss_mixtures_c, ngauss, thisydata, thisycovar, $
      thisamp, thisxmean, thisxcovar, /quiet,$
      _EXTRA=_EXDKEYS
    ;;Collect results
    amps[*,ii]= thisamp
    xmeans[*,*,ii]= thisxmean
    xcovars[*,*,*,ii]= thisxcovar
ENDFOR

;;Now sort the jackknife estimates
IF ndimx EQ 1 THEN BEGIN
    FOR jj=0L, njack -1 DO BEGIN
        sortindx= sort(xmeans[0,*,jj])
        amps[*,jj]= amps[sortindx,jj]
        xmeans[*,*,jj]= xmeans[*,sortindx,jj]
        xcovars[*,*,*,jj]= xcovars[*,*,sortindx,jj]
    ENDFOR
    initsortindx= sort(initxmean)
    amperr= dblarr(ngauss)
    xmeanerr= dblarr(ndimx,ngauss)
    xcovarerr= dblarr(ndimx,ndimx,ngauss)
    FOR kk=0L, ngauss -1 DO BEGIN
        amperr[initsortindx[kk]]= $
          sqrt((njack-1D0)^2D0/njack*variance(amps[kk,*]))
        ;;This only works because dimx=1
        xmeanerr[initsortindx[kk]]= $
          sqrt((njack-1D0)^2D0/njack*variance(xmeans[*,kk,*]))
        xcovarerr[initsortindx[kk]]= $
          sqrt((njack-1D0)^2D0/njack*variance(xcovars[*,*,kk,*]))
    ENDFOR
ENDIF ELSE BEGIN
    print, "ndimx > 1 not implemented at this point"
ENDELSE
END
