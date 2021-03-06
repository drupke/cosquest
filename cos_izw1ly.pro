pro cos_izw1ly

   directoryname = '/Users/drupke/Box Sync/cosquest/fits/'
   gal = 'izw1'
   readcol,directoryname+gal+'/'+gal+'.txt', wavelength, flux, error, $
           FORMAT='D,D,D',/silent

;  Lyalpha
   fittedline = 'Lyalpha'
   contplotreg=[1278.5,1292]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
;   goodind = [[1278.5,1280.7],[1282.3,1283.5],[1284.7,1288.3],[1289.2,1290.2],$
;              [1291.5,1293]]
   goodind = [[1278.5,1280.7],[1282.3,1283.5],[1284.7,1288.3],[1289.2,1289.5],$
              [1289.75,1290.25],[1291,1292]]
;   goodind = [[1278.5,1280.7],[1282.3,1288.3],[1289,1292]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextofit = newind $
      else indextofit = [indextofit,newind]
   endfor
   weight=1d/error^2
;   contfitreg=[[1278.5,1283],[1283,1287],[1287,1288],[1288,1293]]
;   fitfcn=['ifsf_fitpoly','ifsf_fitspline',$
;           'ifsf_fitspline','ifsf_fitspline']
;   fitargs=HASH()
;   fitargs['reg1'] = {}
;   fitargs['reg2'] = {}
;   fitargs['reg3'] = {}
;   fitargs['reg4'] = {}
;   contfitreg=[[1278.5,1285.5],[1285.5,1292]]
;   fitfcn=['ifsf_fitpoly','ifsf_fitspline']
;   fitargs=HASH()
;   fitargs['reg1']={fitord:3}
;   fitargs['reg2']={argsbkpts:{everyn:50}}
;   contfitreg=[[1278.5,1292]]
;   fitfcn=['ifsf_fitspline']
;   fitargs=HASH()
;   fitargs['reg1']={argsbkpts:{everyn:35}}
   contfitreg=[[1278.5,1283],[1283,1286],[1286,1292]]
   fitfcn=['ifsf_fitpoly','ifsf_fitpoly','ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1']={fitord:3}
   fitargs['reg2']={fitord:1}
   fitargs['reg3']={argsbkpts:{everyn:50}}

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
