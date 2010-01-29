;+
;   NAME:
;      leave_one_out
;   PURPOSE:
;      return an lonarr with the i^th  value left out
;   CALLING SEQUENCE:
;      new_array= leave_one_out(array,indx)
;   INPUT:
;      array   - the array
;      indx    - the index of the value you want to leave out
;   OUTPUT:
;      new array
;   REVISION HISTORY:
;      2009-04-06 - Written Bovy (NYU)
;-
FUNCTION LEAVE_ONE_OUT, array, indx

nindx= n_elements(array)
new_array= lonarr(nindx-1)
jj= 0L
FOR ii=0L, nindx-1 DO BEGIN
    IF ii NE indx THEN BEGIN
        new_array[jj]= array[ii]
        jj+= 1
    ENDIF
ENDFOR

RETURN, new_array

END
