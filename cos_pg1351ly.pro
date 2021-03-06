pro cos_pg1351ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1351'
   zsys = 0.0882
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha
;  spline/poly fits
;   fittedline = 'Lyalpha'
;   contplotreg=[1290,1325]
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;                VALUE_LOCATE(wavelength,contplotreg[1])]
;   goodind = [[1290,1300],[1302.7,1303.3],[1305.5,1307],[1311,1312.5],$
;              [1320.5,1323]]
;   for i=0,n_elements(goodind[0,*])-1 do begin
;      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
;                    VALUE_LOCATE(wavelength,goodind[0,i]),$
;                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
;      if i eq 0 then indextofit = newind $
;      else indextofit = [indextofit,newind]
;   endfor
;   weight=1d/error^2
;   contfitreg=[[1290,1312.5],[1311,1322],[1322,1323]]
;   fitfcn=['ifsf_fitspline','ifsf_fitpoly','ifsf_fitpoly']
;   fitargs=HASH()
;   fitargs['reg1'] = {argsbkpts:{everyn:150}}
;   fitargs['reg2'] = {fitord:1}
;   fitargs['reg3'] = {}

;  template fit
   fittedline = 'Lyalpha'
   contplotreg=[1250,1360]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1285,1300],[1305.5,1307],[1311,1312.5],$
              [1320.5,1323]] ;,[1370,1400]]
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
   contfitreg=[[1262,1360]]
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
;      ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
      ifsf_fitmulticont(wavelength, flux, weight, wavelength,itflux, $
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
