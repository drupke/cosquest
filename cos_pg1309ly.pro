pro cos_pg1309ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1309'
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha
   fittedline = 'Lyalpha'
   contplotreg=[1410,1450]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1410,1415],[1416.5,1422],[1424,1427.5],[1428.5,1429],$
              [1429.8,1430.5],[1434.4,1435],[1437,1438.25],[1439.2,1450]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1410,1437.75],[1437.75,1450]]
   fitfcn=['ifsf_fitspline','ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:150}}
   fitargs['reg2'] = {argsbkpts:{everyn:150}}

   set_plot,'z'
   cgplot, wavelength, flux, XRAN=contplotreg, $
           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
           xtit='Wavelength ($\Angstrom$)',$
           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
   continuum=$
      ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
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
   contplotreg=[1204,1213]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1204,1205.6],[1206.8,1207.2],[1207.9,1208.1],[1210.4,1210.6],$
              [1212,1213]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind else indextofit = [indextofit,newind]
   endfor
   weight = 1d/error^2
   contfitreg=[[1204,1213]]
   fitfcn=['ifsf_fitpoly']
   fitargs=HASH()
   fitargs['reg1'] = {fitord:3}

   set_plot,'z'
   cgplot, wavelength, flux, XRAN=contplotreg, $
           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
           xtit='Wavelength ($\Angstrom$)',$
           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
   continuum=$
      ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
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


   ;  Lygamma
   fittedline = 'Lygamma'
   contplotreg=[1145,1151]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
      VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1145,1145.5],[1147.5,1148],[1149.5,1151]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
         VALUE_LOCATE(wavelength,goodind[0,i]),$
         START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind else indextofit = [indextofit,newind]
   endfor
   weight = 1d/error^2
   contfitreg=[[1145,1151]]
   fitfcn=['ifsf_fitpoly']
   fitargs=HASH()
   fitargs['reg1'] = {fitord:1}

   set_plot,'z'
   cgplot, wavelength, flux, XRAN=contplotreg, $
      YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
      1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
      XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
      xtit='Wavelength ($\Angstrom$)',$
      ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
      continuum=$
      ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
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


END
