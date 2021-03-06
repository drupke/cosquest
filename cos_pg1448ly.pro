pro cos_pg1448ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1448'
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha
;  spline/poly
   fittedline = 'Lyalpha'
   contplotreg=[1270,1298]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1270,1277],[1278,1279.5],[1281.5,1282],[1284,1292.3],[1295.5,1298]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1270,1291.5],[1291.5,1298]]
   fitfcn=['ifsf_fitspline' ,'ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:100}}
   fitargs['reg2'] = {argsbkpts:{everyn:100}}

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
