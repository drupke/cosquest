;
; NOTES:
;   width / center unconstrained; definitely uncertain!
;   could possibly constrain center with NV fits
;   

pro cos_pg1001ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1001'
   zsys = 0.161d
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha

; Lorentzian fit
;   fittedline = 'Lyalpha'
;   contplotreg=[1355,1455]
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;                VALUE_LOCATE(wavelength,contplotreg[1])]
;   goodind = [[1355,1360],[1361.5,1365],[1366.5,1370],[1370.8,1371.3],$
;              [1391,1393],[1394.5,1399],[1417,1425],$
;              [1450,1455]]
;   for i=0,n_elements(goodind[0,*])-1 do begin
;      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
;                    VALUE_LOCATE(wavelength,goodind[0,i]),$
;                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
;      if i eq 0 then indextofit = newind $
;      else indextofit = [indextofit,newind]
;   endfor
;   weight=1d/error^2
;   
;   contfitreg=[[1355,1455]]
;   fitfcn=['ifsf_fitpeak']
;   fitargs=HASH()
;   fitargs['reg1'] = {estimates:[2d-14,1411.39,4d,2d-14],$
;                      baseline:1,$
;                      fixed: [0,0,0,0],$
;                      peaktype: 'lorentzian',$
;                      limits: [[0,0],[1405,1415],[1d,10d],[0,0]],$
;                      limited: [[0,0],[1,1],[1,1],[0,0]]}
;; 1215.6701 * 1.161 = 1411.39

;  template fit
   fittedline = 'Lyalpha'
   contplotreg=[1270,1455]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [$
              [1270,1280],[1320,1333],[1336,1350],$
              [1355,1360],[1361.5,1365],[1366.5,1370],[1370.8,1371.3],$
              [1391,1393],[1394.5,1399],[1440,1455]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=dblarr(n_elements(error))
   tolerance = 1d-80
   ibd = where(error^2d LT tolerance)
   igd = where(error^2d GE tolerance)
   weight[ibd] = 0.0d0
   weight[igd] = 1d/error[igd]^2d
   contfitreg=[[1270,1455]]
   fitfcn=['ifsf_fittemplate']
   fitargs=HASH()
   parinfo = replicate({value:0d,fixed:0b},4)
   parinfo[0].value = 1d-14
   fitargs['reg1'] = {templatefcn: 'cos_quasarcomposite',$
                      parinfo:parinfo, npar:4}
;   templatefile = '/Users/drupke/Box Sync/cosquest/spectra/template/'+$
;                  'vandenberk2001.txt'
;   readcol,templatefile,skip=23,/silent,twave,tflux,terr,format='(D,D,D)'
;   templatefile = '/Users/drupke/Box Sync/cosquest/spectra/template/'+$
;                  'agn-composite-stevans-2014.dat'
;   readcol,templatefile,skip=23,/silent,twave,tflux,tspline,format='(D,D,D)'
   templatefile = '/Users/drupke/Box Sync/cosquest/spectra/template/'+$
                  'harris15_composite.fits'
   fxbopen,tlun,templatefile,1
   fxbread,tlun,twave,1
   fxbread,tlun,tflux,2
   fxbclose,tlun
   twave *= 1d + zsys
   itflux = double(ifsf_interptemp(wavelength,twave,tflux))
   itflux /= median(itflux)
   
   set_plot,'z'
   cgplot, wavelength, flux, XRAN=contplotreg, $
;           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;                 1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
           YRAN=[0,3d-14],$
           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
           xtit='Wavelength ($\Angstrom$)',$
           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
   continuum=$
      ifsf_fitmulticont(wavelength, flux, weight, wavelength,itflux, $
                        indextofit,0,fitreg=contfitreg,$
                        fitfcn=fitfcn, fitargs=fitargs,quiet=1b)
   cgoplot, wavelength, continuum, color='Red',thick=4
   img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
                    '_continuum',/jpeg,/nodialog,quality=100)

   nct = contplotind[1]-contplotind[0]+1
   nflux = dblarr(nct)
   nerr = dblarr(nct)
   openw,outlun,directoryname+gal+'/'+gal+fittedline+'_contfit.txt',/get_lun
   for i=0,nct-1 do $
      printf,outlun,wavelength[contplotind[0]+i],flux[contplotind[0]+i],$
             error[contplotind[0]+i],$
             flux[contplotind[0]+i]/continuum[contplotind[0]+i],$
             error[contplotind[0]+i]/continuum[contplotind[0]+i],$
             format='(D12.4,E12.4,E12.4,D8.4,D8.4)'
   free_lun,outlun


;  Lybeta
   fittedline = 'Lybeta'
   
;  poly fit
  contplotreg=[1130,1210]
  contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                 VALUE_LOCATE(wavelength,contplotreg[1])]
goodind = [[1130,1133],[1136,1144],[1146,1152],[1153.5,1160],$
   [1188,1190],[1191,1193],$
   [1194.5,1197.5],[1199,1199.2],[1199.8,1200],$
   [1201.2,1203.5],$
   [1205.7,1206],[1207.7,1210]]
  for i=0,n_elements(goodind[0,*])-1 do begin
     newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                   VALUE_LOCATE(wavelength,goodind[0,i]),$
                   START=VALUE_LOCATE(wavelength,goodind[0,i]))
     if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
  endfor
   weight=dblarr(n_elements(error))
   tolerance = 1d-80
   ibd = where(error^2d LT tolerance)
   igd = where(error^2d GE tolerance)
   weight[ibd] = 0.0d0
   weight[igd] = 1d/error[igd]^2d
;  contfitreg=[[1194.5,1210]]
  contfitreg=[[1130,1157],[1157,1192],[1192,1210]]
  fitfcn=['ifsf_fitspline','ifsf_fitpoly','ifsf_fitspline']
  fitargs=HASH()
  fitargs['reg1'] = {argsbkpts:{everyn:100}}
  fitargs['reg2'] = {}
  fitargs['reg3'] = {argsbkpts:{everyn:100}}
   
   set_plot,'z'
   cgplot, wavelength, flux, XRAN=contplotreg, $
           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
           xtit='Wavelength ($\Angstrom$)',$
           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
   continuum=$
      ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
                        indextoplot,0,fitreg=contfitreg,$
                        fitfcn=fitfcn, fitargs=fitargs,quiet=1b)
   cgoplot, wavelength, continuum, color='Red',thick=4
   img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
                    '_continuum',/jpeg,/nodialog,quality=100)

   
;;  template fit
;   contplotreg=[1070,1220]
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;                VALUE_LOCATE(wavelength,contplotreg[1])]
;   goodind = [[1090,1100],[1120,1130],[1188,1190],[1191,1193],$
;             [1194.5,1197.5],[1199,1199.2],[1199.8,1200],$
;             [1201.2,1203.5],$
;             [1205.7,1206],[1207.7,1210]] ;,[1235,1250]]
;   for i=0,n_elements(goodind[0,*])-1 do begin
;      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
;                    VALUE_LOCATE(wavelength,goodind[0,i]),$
;                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
;      if i eq 0 then indextofit = newind $
;      else indextofit = [indextofit,newind]
;   endfor
;   weight=dblarr(n_elements(error))
;   tolerance = 1d-80
;   ibd = where(error^2d LT tolerance)
;   igd = where(error^2d GE tolerance)
;   weight[ibd] = 0.0d0
;   weight[igd] = 1d/error[igd]^2d
;;  contfitreg=[[1194.5,1210]]
;   contfitreg=[[1070,1220]]
;   fitfcn=['ifsf_fittemplate']
;   fitargs=HASH()
;   parinfo = replicate({value:0d,fixed:0b},4)
;   parinfo[0].value = 1d-14
;   fitargs['reg1'] = {templatefcn: 'cos_quasarcomposite',$
;                      parinfo:parinfo, npar:4}
;   templatefile = '/Users/drupke/Box Sync/cosquest/spectra/template/'+$
;                  'harris15_composite.fits'
;   fxbopen,tlun,templatefile,1
;   fxbread,tlun,twave,1
;   fxbread,tlun,tflux,2
;   fxbclose,tlun
;   twave *= 1d + 0.16d
;   itflux = double(ifsf_interptemp(wavelength,twave,tflux))
;   itflux /= median(itflux)
;   
;   set_plot,'z'
;   cgplot, wavelength, flux, XRAN=contplotreg, $
;;           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;;                 1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
;           YRAN=[0,3d-14],$
;           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
;           xtit='Wavelength ($\Angstrom$)',$
;           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
;   continuum=$
;      ifsf_fitmulticont(wavelength, flux, weight, wavelength,itflux, $
;                        indextofit,0,fitreg=contfitreg,$
;                        fitfcn=fitfcn, fitargs=fitargs,quiet=1b)
;   cgoplot, wavelength, continuum, color='Red',thick=4
;   img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
;                    '_continuum',/jpeg,/nodialog,quality=100)

   nct = contplotind[1]-contplotind[0]+1
   nflux = dblarr(nct)
   nerr = dblarr(nct)
   openw,outlun,directoryname+gal+'/'+gal+fittedline+'_contfit.txt',/get_lun
   for i=0,nct-1 do $
      printf,outlun,wavelength[contplotind[0]+i],flux[contplotind[0]+i],$
             error[contplotind[0]+i],$
             flux[contplotind[0]+i]/continuum[contplotind[0]+i],$
             error[contplotind[0]+i]/continuum[contplotind[0]+i],$
             format='(D12.4,E12.4,E12.4,D8.4,D8.4)'
   free_lun,outlun


   ;  Lygamma
   fittedline = 'Lygamma'
;  poly fit
  contplotreg=[1080,1130]
  contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                 VALUE_LOCATE(wavelength,contplotreg[1])]
  goodind = [[1080,1102],[1115,1130]]
  for i=0,n_elements(goodind[0,*])-1 do begin
     newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                   VALUE_LOCATE(wavelength,goodind[0,i]),$
                   START=VALUE_LOCATE(wavelength,goodind[0,i]))
     if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
  endfor
   weight=dblarr(n_elements(error))
   tolerance = 1d-80
   ibd = where(error^2d LT tolerance)
   igd = where(error^2d GE tolerance)
   weight[ibd] = 0.0d0
   weight[igd] = 1d/error[igd]^2d
;  contfitreg=[[1194.5,1210]]
  contfitreg=[[1080,1099],[1099,1118],[1118,1130]]
  fitfcn=['ifsf_fitpoly','ifsf_fitpoly','ifsf_fitpoly']
  fitargs=HASH()
  fitargs['reg1'] = {}
  fitargs['reg2'] = {fitord:1}
  fitargs['reg3'] = {}
   
   set_plot,'z'
   cgplot, wavelength, flux, XRAN=contplotreg, $
           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
           xtit='Wavelength ($\Angstrom$)',$
           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
   continuum=$
      ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
                        indextoplot,0,fitreg=contfitreg,$
                        fitfcn=fitfcn, fitargs=fitargs,quiet=1b)
   cgoplot, wavelength, continuum, color='Red',thick=4
   img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
                    '_continuum',/jpeg,/nodialog,quality=100)


   nct = contplotind[1]-contplotind[0]+1
   nflux = dblarr(nct)
   nerr = dblarr(nct)
   openw,outlun,directoryname+gal+'/'+gal+fittedline+'_contfit.txt',/get_lun
   for i=0,nct-1 do $
      printf,outlun,wavelength[contplotind[0]+i],flux[contplotind[0]+i],$
             error[contplotind[0]+i],$
             flux[contplotind[0]+i]/continuum[contplotind[0]+i],$
             error[contplotind[0]+i]/continuum[contplotind[0]+i],$
             format='(D12.4,E12.4,E12.4,D8.4,D8.4)'
   free_lun,outlun

END
