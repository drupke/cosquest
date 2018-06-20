pro cos_pg0953ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg0953'
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lybeta
   fittedline = 'Lybeta'
   contplotreg=[1261,1268.5]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1261,1264],[1265.5,1268.5]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1261,1268.5]]
   fitfcn=['ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1']={}

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

;  Lydelta
   fittedline = 'Lydelta'
   contplotreg=[1167,1176]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
      VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1167,1171],[1172.5,1176]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
         VALUE_LOCATE(wavelength,goodind[0,i]),$
         START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1167,1176]]
   fitfcn=['ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1']={}

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
