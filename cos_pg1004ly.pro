pro cos_pg1004ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1004'
   zsys = 0.2406d
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
      FORMAT='D,D,D',/silent

;  Lyalpha
   fittedline = 'Lyalpha'

; Spline fits
;   contplotreg=[1410,1451.8]
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;      VALUE_LOCATE(wavelength,contplotreg[1])]
;   goodind = [[1410,1449],[1451.3,1451.8]]
;   for i=0,n_elements(goodind[0,*])-1 do begin
;      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
;         VALUE_LOCATE(wavelength,goodind[0,i]),$
;         START=VALUE_LOCATE(wavelength,goodind[0,i]))
;      if i eq 0 then indextofit = newind $
;      else indextofit = [indextofit,newind]
;   endfor
;   weight=1d/error^2
;   contfitreg=[[1410,1451.8]]
;   fitfcn=['ifsf_fitspline']
;   fitargs=HASH()
;   fitargs['reg1'] = {argsbkpts:{everyn:100}}
;   set_plot,'z'
;   cgplot, wavelength, flux, XRAN=contplotreg, $
;      YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;      1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
;      XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
;      xtit='Wavelength ($\Angstrom$)',$
;      ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
;      continuum=$
;      ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
;      indextofit,0,fitreg=contfitreg,$
;      fitfcn=fitfcn, fitargs=fitargs)
;   cgoplot, wavelength, continuum, color='Red',thick=4
;   img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
;      '_continuum',/jpeg,/nodialog,quality=100)

; Template fit
   contplotreg=[1310,1450]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
      VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1310,1330],[1340,1390],[1405,1410]]
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
   contfitreg=[[1310,1450]]
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

;  Lybeta
   fittedline = 'Lybeta'
   contplotreg=[1155,1450]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
      VALUE_LOCATE(wavelength,contplotreg[1])]
  goodind = [[1155,1165],[1172,1180],[1228,1232],$
             [1280.2,1280.4],[1281.0,1281.8],[1284.0,1284.2],[1286.3,1286.6],$
             [1287.3,1300],[1308,1332],[1338,1390],[1405,1410]]
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
   contfitreg=[[1155,1450]]
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
;     YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;     1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
     YRAN=[0,2d-14],$
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


;  Lygamma
   fittedline = 'Lygamma'
   contplotreg=[1155,1450]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
      VALUE_LOCATE(wavelength,contplotreg[1])]
  goodind = [[1155,1165],[1172,1180],[1228,1232],$
             [1280.2,1280.4],[1281.0,1281.8],[1284.0,1284.2],[1286.3,1286.6],$
             [1287.3,1300],[1308,1332],[1338,1390],[1405,1410]]
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
   contfitreg=[[1155,1450]]
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
;     YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;     1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
     YRAN=[0,2d-14],$
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
