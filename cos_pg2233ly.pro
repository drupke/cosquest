pro cos_pg2233ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg2233'
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lybeta
;  spline/poly
   fittedline = 'Lybeta'
   contplotreg=[1350,1380]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
      VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1350,1354],[1355.5,1357],[1358.5,1359],[1360,1360.5],$
      [1363.5,1366],[1371,1373],[1374.5,1380]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
         VALUE_LOCATE(wavelength,goodind[0,i]),$
         START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1350,1380]]
   fitfcn=strarr(n_elements(contfitreg[0,*]))+'ifsf_fitspline'
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:200}}

   set_plot,'z'
   cgplot, wavelength, flux, XRAN=contplotreg, $
           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
           1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
           xtit='Wavelength ($\Angstrom$)',$
           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
   continuum=ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
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
