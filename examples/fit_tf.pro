;+
;   NAME:
;      fit_TF
;   PURPOSE:
;      fit the Tully-Fisher relation from HST cepheids using using 1
;      Gaussian mixture
;   CALLING SEQUENCE:
;   INPUT:
;      plotfilename  - filename for plot
;      texfilename   - filename for results of fit
;   OUTPUT:
;   REVISION HISTORY:
;      2009-04-06 - Written Bovy (NYU)
;-
PRO FIT_TF, plotfilename=plotfilename, texfilename=texfilename

IF ~keyword_set(plotfilename) THEN plotfilename='TF.ps'
IF ~keyword_set(texfilename) THEN texfilename='TF.tex'
tol= 1D-5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; READ AND SET UP THE DATA
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;Read the data
missing_value= 10000000.0
nfield= 19L
; most data types are double
fieldtypes= lonarr(nfield)+5
; some are long
fieldtypes[[15,16,17,18]]= 3
; some are string
fieldtypes[[0]]= 7
; give dumb names
fieldnames= 'TF'+string(indgen(nfield),format='(I2.2)')
template= {version:1.0, datastart: 0, $
           delimiter: '|', missingvalue: missing_value, commentsymbol: '#', $
           fieldcount: nfield, fieldtypes: fieldtypes, $
           fieldnames: fieldnames, $
           fieldlocations: lonarr(nfield), fieldgroups: lindgen(nfield)}
tfdata= read_ascii('TF.dat',num_records=nline,template=template)
print, 'found '+strtrim(string(n_elements(tfdata.tf01)),2)+$
  ' records'

;;Create the data structure
datastruct={TF, $
            name     :' '       , $  ; Name of the galaxy
            mags     : dblarr(5), $  ; magnitudes in different bands [B,V,R,I,H]
            mags_err : dblarr(5), $  ; errors in magnitudes
            W20      : 0D       , $  ; log of W (20%)
            W20_err  : 0D       , $  ; error in W20
            W50      : 0D       , $  ; log of W (50%)
            W50_err  : 0D       , $  ; error in W50
            incI     : 0L       , $  ; inclination
            incI_err : 0L       , $  ; error iI
            inc      : 0L       , $  ; inclination
            inc_err  : 0L         $  ; error in i
}
tf= replicate(datastruct,n_elements(tfdata.tf01))
;;load the data in the structure
tf.name=          tfdata.tf00
tf.mags[0]=       tfdata.tf01
tf.mags[1]=       tfdata.tf03
tf.mags[2]=       tfdata.tf05
tf.mags[3]=       tfdata.tf07
tf.mags[4]=       tfdata.tf09
tf.mags_err[0]=   tfdata.tf02
tf.mags_err[1]=   tfdata.tf04
tf.mags_err[2]=   tfdata.tf06
tf.mags_err[3]=   tfdata.tf08
tf.mags_err[4]=   tfdata.tf10
tf.W20=           tfdata.tf11
tf.W20_err=       tfdata.tf12
tf.W50=           tfdata.tf13
tf.W50_err=       tfdata.tf14
tf.incI=          tfdata.tf15
tf.incI_err=      tfdata.tf16
tf.inc=           tfdata.tf17
tf.inc_err=       tfdata.tf18


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; MAIN FIT OF TF RELATION
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
slopes= dblarr(5)
slopes_err= dblarr(5)
intercepts= dblarr(5)
intercepts_err= dblarr(5)
FOR ii=0L, 4 DO BEGIN
    ;;Set-up the various arrays
    ;;First check whether any magnitudes are missing
    indx= where(tf.mags[ii] NE missing_value)
    ndata= n_elements(indx)
    ydata= dblarr(2,ndata)
    ycovar= dblarr(2,2,ndata)
    projection= dblarr(2,2,ndata)
    amp= dblarr(1)+1.
    xmean= dblarr(2,1)
    xcovar= dblarr(2,2,1)
    ydata[0,0:ndata-1]= tf[indx].W20
    ydata[1,0:ndata-1]= tf[indx].mags[ii]
    ycovar[0,0,0:ndata-1]= (tf[indx].W20_err)^2.
    ycovar[1,1,0:ndata-1]= (tf[indx].mags_err[ii])^2.
    FOR jj=0L, ndata-1 DO BEGIN
        projection[0,0,jj]= 1.
        projection[1,1,jj]= 1.
    ENDFOR
    xcovar[0,0,0]= (max(ydata[0,0:ndata-1])-min(ydata[0,0:ndata-1]))^2.
    xcovar[1,1,0]= (max(ydata[1,0:ndata-1])-min(ydata[1,0:ndata-1]))^2.
    ;;Run proj_gauss_mixtures
    projected_gauss_mixtures_c, 1, ydata, ycovar, $
      amp, xmean, xcovar, tol=tol, /quiet, projection=projection
    ;;Slope and zero point
    eigenvals= EIGENQL(xcovar,eigenvectors=eigenvectors)
    slopes[ii]= eigenvectors[1,0]/eigenvectors[0,0]
    intercepts[ii]= 2.5*slopes[ii]-slopes[ii]*xmean[0]+xmean[1]
    ;;Jackknife the values of the slope and intercept
    jack_slopes= dblarr(ndata)
    jack_int= dblarr(ndata)
    FOR jj= 0L, ndata-1 DO BEGIN
        ;;Create new sample by leaving one datapoint out
        jackindx= leave_one_out(indx,jj)
        jack_ydata= dblarr(2,ndata-1)
        jack_ycovar= dblarr(2,2,ndata-1)
        jack_projection= dblarr(2,2,ndata-1)
        amp= dblarr(1)+1.
        xmean= dblarr(2,1)
        xcovar= dblarr(2,2,1)
        jack_ydata[0,0:ndata-2]= tf[jackindx].W20
        jack_ydata[1,0:ndata-2]= tf[jackindx].mags[ii]
        jack_ycovar[0,0,0:ndata-2]= (tf[jackindx].W20_err)^2.
        jack_ycovar[1,1,0:ndata-2]= (tf[jackindx].mags_err[ii])^2.
        FOR kk=0L, ndata-2 DO BEGIN
            jack_projection[0,0,kk]= 1.
            jack_projection[1,1,kk]= 1.
        ENDFOR
        xcovar[0,0,0]= (max(ydata[0,0:ndata-2])-min(ydata[0,0:ndata-2]))^2.
        xcovar[1,1,0]= (max(ydata[1,0:ndata-2])-min(ydata[1,0:ndata-2]))^2.
        ;;Run proj_gauss_mixtures
        projected_gauss_mixtures_c, 1, jack_ydata, jack_ycovar, $
          amp, xmean, xcovar, tol=tol, /quiet, projection=jack_projection
        ;;Slope and zero point
        eigenvals= EIGENQL(xcovar,eigenvectors=eigenvectors)
        jack_slopes[jj]= eigenvectors[1,0]/eigenvectors[0,0]
        jack_int[jj]= 2.5*jack_slopes[jj]-jack_slopes[jj]*xmean[0]+xmean[1]
    ENDFOR
    squared_err= (ndata-1.)^2./double(ndata)*VARIANCE(jack_slopes,/double)
    slopes_err[ii]= sqrt(squared_err)
    squared_err= (ndata-1.)^2./double(ndata)*VARIANCE(jack_int,/double)
    intercepts_err[ii]= sqrt(squared_err)
ENDFOR

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; WRITE THE RESULTS TO A FILE
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
bandnames= ['B_T^c','V_T^c','R^c_T','I^c_T','H^c_{-0.5}']
xname='\log_{10} W_{20}^c'
OPENW, wlun, texfilename, /GET_LUN
PRINTF, wlun, '\begin{align}'
FOR ii=0L, 3 DO BEGIN
    PRINTF, wlun, bandnames[ii]+' &= -('+$
      strtrim(string(abs(slopes[ii]),format='(F4.2)'),2)+$
      ' \pm '+strtrim(string(abs(slopes_err[ii]),format='(F4.2)'),2)+')\,( '+$
      xname+' - 2.5 ) - ( '+$
      strtrim(string(abs(intercepts[ii]),format='(F5.2)'),2)+$
      ' \pm '+strtrim(string(abs(intercepts_err[ii]),format='(F4.2)'),2)+$
      ') \\'
ENDFOR
ii=4
PRINTF, wlun, bandnames[ii]+' &= -('+$
  strtrim(string(abs(slopes[ii]),format='(F5.2)'),2)+$
  ' \pm '+strtrim(string(abs(slopes_err[ii]),format='(F4.2)'),2)+')\,( '+$
  xname+' - 2.5 ) - ( '+$
  strtrim(string(abs(intercepts[ii]),format='(F5.2)'),2)+$
  ' \pm '+strtrim(string(abs(intercepts_err[ii]),format='(F4.2)'),2)+$
  ') '
PRINTF, wlun, '\end{align}'
FREE_LUN, wlun


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; PLOT THE RESULT
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;Plot everything
bands=['B','V','R','I','H']
xrange=[2.25,2.85]
yranges= dblarr(2,5)
yranges[0:1,0]=[-17.5,-23]
yranges[0:1,1]=[-18,-23.5]
yranges[0:1,2]=[-18.3,-23.75]
yranges[0:1,3]=[-18.75,-24.22]
yranges[0:1,4]=[-19.55,-25]
xtitle='log!d10!n W (20%)'

; setup postscript file
if(NOT keyword_set(axis_char_scale)) then axis_char_scale= 1.75
if(NOT keyword_set(tiny)) then tiny=1.d-4
pold=!P
xold=!X
yold=!Y
bangp=!P
bangx=!X
bangy=!Y
!P.FONT= -1
set_plot, "PS"
!P.BACKGROUND= 16777215
!P.COLOR= 0
if(NOT keyword_set(xsize)) then xsize= 3.375*2.
if(NOT keyword_set(ysize)) then ysize= 6.
device, file=plotfilename,/inches,xsize=xsize,ysize=ysize, $
  xoffset=(8.5-xsize)/2.0,yoffset=(11.0-ysize)/2.0,/color, $
  bits_per_pixel=64
!P.THICK= 1.0
!P.CHARTHICK= !P.THICK & !X.THICK= !P.THICK & !Y.THICK= !P.THICK
!P.CHARSIZE= 0.5
!P.PSYM= 0
!P.LINESTYLE= 0
!P.TITLE= ''
!X.STYLE= 1
!X.CHARSIZE= axis_char_scale
!X.MARGIN= [1,1]*0.5
!X.OMARGIN= [7,7]*axis_char_scale
!X.RANGE= 0
!X.TICKS= 0
!Y.STYLE= 1
!Y.CHARSIZE= !X.CHARSIZE
!Y.MARGIN= 0.6*!X.MARGIN
!Y.OMARGIN= 0.6*!X.OMARGIN
!Y.RANGE= 0
!Y.TICKS= !X.TICKS
!P.MULTI= [1,1,1]
xyouts, 0,0,'!6'
colorname= ['red','green','blue','magenta','cyan','dark yellow', $
            'purple','light green','orange','navy','light magenta', $
            'yellow green']
ncolor= n_elements(colorname)
loadct,0



;;plotting symbol
phi=findgen(32)*(!PI*2/32.)
phi = [ phi, phi(0) ]
usersym, cos(phi), sin(phi), /fill
xpos_label= 2.35
ypos_label_rel= (yranges[0,0]+22.4)/(yranges[0,0]-yranges[1,0])

positions= dblarr(4,5)
positions[0:3,0]= [0.08,0.55,0.306,0.9] 
positions[0:3,1]= [0.386,0.55,0.616,0.9]
positions[0:3,2]= [0.696,0.55,0.924,0.9] 
positions[0:3,3]= [0.234,0.1,0.46,0.45] 
positions[0:3,4]= [0.54,0.1,0.76,0.45] 

FOR ii=0L, 4 DO BEGIN
    xline= [xrange[0],xrange[1]]
    yline= [slopes[ii]*(xline[0]-2.5)+intercepts[ii],slopes[ii]*(xline[1]-2.5)+intercepts[ii]]
    ploterror, tf.W20, tf.mags[ii], tf.W20_err, tf.mags_err[ii], psym=8, yrange=yranges[0:1,ii], $
      xrange=xrange, xtitle=xtitle, position=positions[0:3,ii], symsize=0.5, /NOERASE
    oplot, xline, yline, psym=-3
    ypos_label=yranges[0,ii]+ypos_label_rel*(yranges[1,ii]-yranges[0,ii])
    hasLegend= 0
    CATCH, Error_status
    IF Error_status NE 0 THEN BEGIN  
        print, "Skipping legend because legend function not found"
        hasLegend= 1
        CATCH, /CANCEL  
    ENDIF  
    IF hasLegend EQ 0 THEN legend, [bands[ii]],pos=[xpos_label,ypos_label], box=0, charsize=1.5*!P.charsize
ENDFOR

device,/close


END
