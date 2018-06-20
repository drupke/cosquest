pro cos_pg1116ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1116'
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha
   fittedline = 'Lyalpha'
   contplotreg=[1380,1430]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1380,1383],[1385,1393],[1395,1402],$
              [1404,1416],[1419.5,1425.5],[1427.5,1430]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1380,1430]]
   fitfcn=['ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:100}}

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


;  Lybeta
   fittedline = 'Lybeta'
   contplotreg=[1165,1214]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1165,1167.5],[1168.5,1174.5],[1175.5,1179.5],[1180.5,1181],$
              [1182,1189.5],[1191.5,1192.5],[1194.5,1195.5],$
              [1197.5,1198.7],[1201.5,1202],[1203,1203.2],$
              [1204.3,1205.7],[1207.8,1208.8],[1210.2,1210.5],[1211.2,1214]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1165,1214]]
   fitfcn=['ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:100}}

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
