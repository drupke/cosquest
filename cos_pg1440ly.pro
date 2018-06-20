pro cos_pg1440ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1440'
   zsys = 0.077d
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha
;  spline/poly
   fittedline = 'Lyalpha'
   contplotreg=[1290,1315]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1290,1294.5],[1295.5,1299],[1300.5,1301],[1302.5,1303.5],$
              [1304.5,1305.5],[1306,1306.5],[1307.5,1309.5],[1310.5,1315]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1290,1298],[1298,1315]] ;,[1308,1315]]
   fitfcn=['ifsf_fitpoly','ifsf_fitspline' ,'ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1'] = {}
   fitargs['reg2'] = {argsbkpts:{everyn:100}}
;   fitargs['reg3'] = {argsbkpts:{everyn:100}}

; template
;   fittedline = 'Lyalpha'
;   contplotreg=[1275,1320]
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;                VALUE_LOCATE(wavelength,contplotreg[1])]
;   goodind = [[1275,1294.5],[1295.5,1299],[1300.5,1301],[1302.5,1303.5],$
;              [1304.5,1305.5],[1306,1306.5],[1307.5,1309.5],[1310.5,1320]]
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
;   contfitreg=[[1275,1320]]
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
;   twave *= 1d + zsys
;   itflux = double(ifsf_interptemp(wavelength,twave,tflux))
;   itflux /= median(itflux)

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
