;+
;   NAME:
;      weighted_kstwo
;   PURPOSE:
;      perform a weighted KS test
;   CALLING SEQUENCE:
;      weighted_kstwo, data1, data2, weight1=weight1, weight2=weight2, $
;      D, prob
;   INPUT:
;       data1 -  vector of data values, at least 4 data values must be
;               included for the K-S statistic to be meaningful
;       data2 -  second set of data values, does not need to have the
;               same number of elements as data1
;   OPTIONAL INPUT:
;      weight1 - weights for the first data vector
;      weight2 - weights for the second data vector
;   OUTPUT:
;      D - floating scalar giving the Kolmogorov-Smirnov statistic.
;           It specifies the maximum deviation between the cumulative 
;            distribution of the data and the supplied function 
;      prob - floating scalar between 0 and 1 giving the significance
;             level of the K-S statistic.   Small values of PROB show that
;             the cumulative distribution function of DATA1 is significantly 
;             different from DATA2
;   HISTORY:
;      2010-03-22 - Written (adapted from kstwo) - Bovy (NYU)
;-
PRO WEIGHTED_KSTWO, data1, data2, weight1=weight1, weight2=weight2, $
                    D, prob
On_error, 2
compile_opt idl2 
  
if ( N_params() LT 4 ) then begin
    print,'Syntax - WEIGHTED_KSTWO, data1, data2, d, prob, weighted1=weighted1, weighted2=weighted2'
    return
endif

n1 = N_elements( data1 )
if ( n1 LE 3 ) then message, $
  'ERROR - Input data values (first param) must contain at least 4 values'

n2 = N_elements( data2 )
if ( n2 LE 3 ) then message, $
  'ERROR - Input data values (second param) must contain at least 4 values'

sortindx1= sort(data1)
sortdata1 = data1[ sortindx1 ] ;Sort input arrays into 
sortindx2= sort(data2)
sortdata2 = data2[ sortindx2 ] ;ascending order

IF ~keyword_set(weight1) THEN BEGIN
    fn1 = ( dindgen( n1 +1 )  ) / n1
ENDIF ELSE BEGIN
    fn1= total(weight1[sortindx1],/cumulative)/total(weight1)
ENDELSE
IF ~keyword_set(weight2) THEN BEGIN
    fn2 = ( dindgen( n2 +1)  ) / n2
ENDIF ELSE BEGIN
    fn2= total(weight2[sortindx2],/cumulative)/total(weight2)
ENDELSE

j1 = 0l & j2 = 0l
id1 = lonarr(n1+n2)  & id2 = id1
i = 0l

; Form the two cumulative distribution functions, marking points where
; one
; must test their difference

while ( j1 LT n1 ) and ( j2 LT n2 ) do begin
    
    d1 = sortdata1[j1]
    d2 = sortdata2[j2]
    if d1 LE d2 then j1 = j1 +1
    if d2 LE d1 then j2 = j2 +1
    
    id1[i] = j1   & id2[i] = j2
    i = i+1
    
endwhile

id1 = id1[0:i-1]   &  id2 = id2[0:i-1]

; The K-S statistic D is the maximum difference between the two
; distribution
; funtions

D = max( abs( fn1[id1] - fn2[id2] ) ) 
N_eff =  n1*n2/ float(n1 + n2)  ;Effective # of data points
PROB_KS, D, N_eff, prob         ;Compute significance of statistic

return
end
