pro cos_pg1617ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1617'
   zsys=0.114d
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha
   fittedline = 'Lyalpha'
  contplotreg=[1280,1440]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
      VALUE_LOCATE(wavelength,contplotreg[1])]
  goodind = [[1141,1142],[1147.5,1151.5],[1153,1189],[1191,1192.5],[1194,1198],$
   [1201.5,1205],[1207.5,1224.5],[1225.5,1231.5],[1232.5,1236.5],$
   [1237.7,1238.2],[1239.2,1239.3],[1240.5,1242],[1243.5,1250],$
   [1251.5,1253.5],[1254.5,1257.5],[1261,1276.5],[1278.5,1301],$
   [1305,1320],[1321,1325],$
   [1326,1327],[1330,1333.5],[1336.5,1339],$
   [1344,1346.5],[1347.5,1349],[1350.5,1360],[1369,1369.8],[1374,1392.5],$
   [1394.5,1402],$
   [1404,1440]]
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
  contfitreg=[[1280,1440]]
   fitfcn=['ifsf_fittemplate']
   fitargs=HASH()
   parinfo = replicate({value:0d,fixed:0b},4)
   parinfo[0].value = 1d-14
   fitargs['reg1'] = {templatefcn: 'cos_quasarcomposite',$
      parinfo:parinfo, npar:4}
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
           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
      XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
      xtit='Wavelength ($\Angstrom$)',$
      ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
      continuum=$
      ifsf_fitmulticont(wavelength, flux, weight, wavelength, itflux, $
      indextofit,0,fitreg=contfitreg,$
      fitfcn=fitfcn, fitargs=fitargs)
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


;  spline/poly
;   fittedline = 'Lyalpha'
;   contplotreg=[1330,1354]
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;                VALUE_LOCATE(wavelength,contplotreg[1])]
;   goodind = [[1330,1333.5],[1336.5,1339],[1344,1346.5],[1347.5,1349],$
;      [1350.5,1354]]
;   for i=0,n_elements(goodind[0,*])-1 do begin
;      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
;                    VALUE_LOCATE(wavelength,goodind[0,i]),$
;                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
;      if i eq 0 then indextofit = newind $
;      else indextofit = [indextofit,newind]
;   endfor
;   weight=1d/error^2
;   contfitreg=[1330,1354]
;   fitfcn=['ifsf_fitpoly']
;   fitargs=HASH()
;   fitargs['reg1'] = {}
;
;   set_plot,'z'
;   cgplot, wavelength, flux, XRAN=contplotreg, $
;           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
;           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
;           xtit='Wavelength ($\Angstrom$)',$
;           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
;   continuum=$
;      ifsf_fitmulticont(wavelength, flux, weight, wavelength, itflux, $
;                        indextofit,0,fitreg=contfitreg,$
;                        fitfcn=fitfcn, fitargs=fitargs)
;   cgoplot, wavelength, continuum, color='Red',thick=4
;   img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
;                    '_continuum',/jpeg,/nodialog,quality=100)
;
;   nct = contplotind[1]-contplotind[0]+1
;   nflux = dblarr(nct)
;   nerr = dblarr(nct)
;   openw,outlun,directoryname+gal+'/'+gal+fittedline+'_contfit.txt',/get_lun
;   for i=0,nct-1 do $
;      printf,outlun,wavelength[contplotind[0]+i],flux[contplotind[0]+i],$
;             error[contplotind[0]+i],$
;             flux[contplotind[0]+i]/continuum[contplotind[0]+i],$
;             error[contplotind[0]+i]/continuum[contplotind[0]+i],$
;             format='(D12.4,E12.4,E12.4,D8.4,D8.4)'
;   free_lun,outlun

   ;  Lybeta
   fittedline = 'Lybeta'
  contplotreg=[1133,1175]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
      VALUE_LOCATE(wavelength,contplotreg[1])]
  goodind = [[1141,1142],[1147.5,1151.5],[1153,1189],[1191,1192.5],[1194,1198],$
   [1201.5,1205],[1207.5,1224.5],[1225.5,1231.5],[1232.5,1236.5],$
   [1237.7,1238.2],[1239.2,1239.3],[1240.5,1242],[1243.5,1250],$
   [1251.5,1253.5],[1254.5,1257.5],[1261,1276.5],[1278.5,1301],$
   [1305,1320],[1321,1325],[1326,1327],[1330,1333.5],[1336.5,1339],$
   [1344,1346.5],[1347.5,1349],[1350.5,1360],[1369,1369.8],[1374,1390]]
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
  contfitreg=[[1133,1175]]
   fitfcn=['ifsf_fittemplate']
   fitargs=HASH()
   parinfo = replicate({value:0d,fixed:0b},4)
   parinfo[0].value = 1d-14
   fitargs['reg1'] = {templatefcn: 'cos_quasarcomposite',$
      parinfo:parinfo, npar:4}
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
           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
      XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
      xtit='Wavelength ($\Angstrom$)',$
      ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
      continuum=$
      ifsf_fitmulticont(wavelength, flux, weight, wavelength, itflux, $
      indextofit,0,fitreg=contfitreg,$
      fitfcn=fitfcn, fitargs=fitargs)
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
