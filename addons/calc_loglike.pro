;+
;   NAME:
;      calc_loglike
;   PURPOSE:
;      calculate the likelihood of a given point given a decompostion
;      of a probability distribution as a sum of Gaussians
;   CALLING SEQUENCE:
;      loglike= calc_loglike(ydata, ycovar, xmean, xcovar, xamp, projection=projection)
;   INPUT:
;      ydata    - the data [ndimy,ndata]
;      ycovar   - error covariances [ndimy,ndimy,ndata] or
;                 [ndimy,ndata] for diagonal covariances
;      xmean     - means of the gaussians [ndimx,ngauss]
;      xcovar    - covariances of the gaussians [ndimx,ndimx,ngauss]
;      xamp      - amplitudes of the gaussians [ngauss]
;   OPTIONAL INPUTS:
;      projection - [ndimx, ndimy, ndata] non-square matrices
;                  implementing a projection of model space onto data
;                  space (useful if your data occupy a
;                  lower-dimensional subspace of the model space, or
;                  live in a rotated version of model space, ...)
;   OUTPUT:
;      loglike     - value or array of ln of posterior probabilities [ndata]
;   REVISION HISTORY:
;      2008-12-16 - Written - Jo Bovy (NYU)
;      2010-02-24 - Re-written to be stand-alone - Bovy
;      2010-03-10 - Calculate log-likelihood, not posterior
;                   probability to belong to one of the Gaussians - Bovy
;-
FUNCTION BOVY_LOGSUM, array
amax= max(array)
RETURN, amax+alog(total(exp(array-amax),/double))
END
FUNCTION SPECIAL_INVERT, matrix
RETURN, invert(matrix,/double)
END
FUNCTION BOVY_DETERM, matrix, double=double, check=check
IF n_elements(matrix) EQ 1 THEN return, matrix ELSE $
  return, determ(matrix,double=double,check=check)
END
FUNCTION CALC_LOGLIKE, ydata, ycovar, xmean, xcovar, xamp, $
                       projection=projection

ngauss= n_elements(xamp)
ndimx= n_elements(xmean)/ngauss
ndimy=(size(ydata,/dimensions))[0]
ndata=n_elements(ydata)/ndimy
if ~keyword_set(projection) and ndimx NE ndimy THEN $
  message, "Dimension of model is not equal to the dimension of the data and projection is not set"
IF n_elements(ycovar) EQ n_elements(ydata) THEN diagcovar= 1B ELSE $
  diagcovar= 0B

twopiterm=0.5*double(ndimy)*alog(2.*!DPI)

IF ndata EQ 1 THEN scalarout= 1B ELSE scalarout= 0B

out= dblarr(ndata)

;;Loop over data and Gaussians to find posterior probabilities
FOR ii= 0L, ndata-1 DO BEGIN
    loglike= dblarr(ngauss)
    FOR kk= 0L, ngauss-1 DO BEGIN
        IF keyword_set(projection) THEN BEGIN
            IF diagcovar THEN BEGIN
                tinv=special_invert(transpose(projection[*,*,ii])#xcovar[*,*,kk]# $
                                  projection[*,*,ii]+diag_matrix(ycovar[*,ii]))
            ENDIF ELSE BEGIN
                tinv=special_invert(transpose(projection[*,*,ii])#xcovar[*,*,kk]# $
                                    projection[*,*,ii]+ycovar[*,*,ii])
            ENDELSE
            delta=ydata[*,ii]-transpose(projection[*,*,ii])#xmean[*,kk]
        ENDIF ELSE BEGIN
            IF diagcovar THEN BEGIN
                tinv=special_invert(xcovar[*,*,kk]+diag_matrix(ycovar[*,ii]))
            ENDIF ELSE BEGIN
                tinv=special_invert(xcovar[*,*,kk]+ycovar[*,*,ii])
            ENDELSE
            delta=ydata[*,ii]-xmean[*,kk]
        ENDELSE
        loglike[kk]=alog(xamp[kk])+0.5*alog(bovy_determ(tinv,/double,/check) > (machar(/double)).xmin)- $
          0.5*transpose(delta)#tinv#delta-twopiterm
    ENDFOR
    ;;sum the probabilities under the invidual components
    out[ii]= bovy_logsum(loglike)
ENDFOR
IF scalarout THEN RETURN, out[0] ELSE RETURN, out
END
