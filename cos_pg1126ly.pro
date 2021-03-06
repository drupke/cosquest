pro cos_pg1126ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'pg1126'
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha
   fittedline = 'Lyalpha'
   contplotreg=[1265,1293]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1265,1267],[1267.5,1269],[1272.5,1275],[1277.75,1277.95],[1282,1283],$
              [1285.5,1286.2],[1287.1,1287.8],[1290,1293]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
   contfitreg=[[1265,1287.45],[1287.45,1293]]
   fitfcn=['ifsf_fitspline','ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:100}}
   fitargs['reg2'] = {argsbkpts:{everyn:50}}

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
   
   readcol,directoryname+gal+'/'+gal+'OVI.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

   fittedline = 'Lybeta'
   contplotreg=[1066,1120]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   goodind = [[1066,1068.5],[1074,1075],[1080,1082],[1088.4,1088.6],$
              [1094.4,1095.1],[1097.5,1098],[1098.6,1099],[1101,1107.5],$
              [1109,1109.5],[1110.5,1111.5],[1113,1115],[1115.5,1120]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind else indextofit = [indextofit,newind]
   endfor
   weight = 1d/error^2
   contfitreg=[[1066,1075],[1074,1095.8],[1094.8,1097.75],[1097.75,1120]]
   fitfcn=['ifsf_fitpoly','ifsf_fitpoly','ifsf_fitpoly','ifsf_fitpoly']
   fitargs=HASH()
   fitargs['reg1'] = {fitord:1}
   fitargs['reg2'] = {fitord:2}
   fitargs['reg3'] = {fitord:2}
   fitargs['reg4'] = {fitord:20}

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
