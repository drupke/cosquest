FUNCTION cos_pg1004ovi, directoryname, gal, zgal, profileshifts, $
                        profilesig, coveringfactor,opticaldepth

  fittedline = 'OVI'
  bad = 1d99
  gal = gal[0]
  bin = 2d
  ncols = 1
  nrows = 1
  centcol = 1
  centrow = 1
  zgal=zgal[0]
  outstr = 'rb'+string(bin,format='(I0)')
  comps=N_ELEMENTS(profileshifts)
  readcol,directoryname+'/'+gal+'/'+gal+'.txt', wavelength, flux, error, $
          FORMAT='D,D,D',/silent
  
;  Fitting low-v OVI only
;; Finding the index to fit over
;   linefitreg=[1265,1288]
;   lineplotreg=[1263,1292]
;   contplotreg=[1249,1295]
;   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;                VALUE_LOCATE(wavelength,contplotreg[1])]
;   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
;               VALUE_LOCATE(wavelength,linefitreg[1])]
;   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
;                VALUE_LOCATE(wavelength,lineplotreg[1])]
;  goodind = [[1249,1250.2],[1251,1252],[1252.8,1253.3],[1254.5,1259],$
;             [1262.5,1265.5],[1266.5,1267],[1268,1270],$
;             [1271,1271.9],[1273.5,1274.5],[1279.4,1279.7],[1280.2,1280.4],$
;             [1281.0,1281.8],[1284.0,1284.2],[1286.3,1286.6],[1287.3,1295]]
;  for i=0,n_elements(goodind[0,*])-1 do begin
;     newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
;        VALUE_LOCATE(wavelength,goodind[0,i]),$
;        START=VALUE_LOCATE(wavelength,goodind[0,i]))
;     if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
;  endfor
;  weight=1d/error^2
;  contfitreg=[[1249,1256],[1256,1279.55],[1279.55,1286.55],[1286.55,1295]]
;  fitfcn=['ifsf_fitspline','ifsf_fitpoly','ifsf_fitpoly','ifsf_fitpoly']
;  fitargs=HASH()
;  fitargs['reg1'] = {argsbkpts:{everyn:100}}
;  fitargs['reg2'] = {fitord:2}
;  fitargs['reg3'] = {fitord:3}
;  fitargs['reg4'] = {fitord:2}
;;  fitargs['reg3'] = {argsbkpts:{everyn:100}}

; Finding the index to fit over
   linefitreg=[1235,1288]
   lineplotreg=[1230,1292]
   contplotreg=[1155,1450]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
               VALUE_LOCATE(wavelength,linefitreg[1])]
   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
                VALUE_LOCATE(wavelength,lineplotreg[1])]
  goodind = [[1155,1165],[1172,1180],[1228,1232],$
             [1280.2,1280.4],[1281.0,1281.8],[1284.0,1284.2],[1286.3,1286.6],$
             [1287.3,1300],[1308,1332],[1338,1390],[1405,1410]]
  for i=0,n_elements(goodind[0,*])-1 do begin
     newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
        VALUE_LOCATE(wavelength,goodind[0,i]),$
        START=VALUE_LOCATE(wavelength,goodind[0,i]))
     if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
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
   twave *= 1d + zgal[0]
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
     continuum=ifsf_fitmulticont(wavelength, flux, weight, wavelength,itflux, $
     indextoplot,0,fitreg=contfitreg,$
     fitfcn=fitfcn, fitargs=fitargs)
  cgoplot, wavelength, continuum, color='Red',thick=4
  img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
     '_continuum',/jpeg,/nodialog,quality=100)

  relativeflux=flux/continuum
  relativeerror=error/continuum
  
  ; Removing problem areas by setting error to 0 so that MPFITFUN ignores it
  relativeerror[VALUE_LOCATE(wavelength,1250):$
                VALUE_LOCATE(wavelength,1251)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1251.7):$
                VALUE_LOCATE(wavelength,1253)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1253.5):$
                VALUE_LOCATE(wavelength,1255)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1259):$
                VALUE_LOCATE(wavelength,1262.5)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1267):$
                VALUE_LOCATE(wavelength,1272.5)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1280.3):$
                VALUE_LOCATE(wavelength,1281.2)] = bad
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Parameters for doublet fit
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ndoubletabs = dblarr(ncols,nrows)+comps
    doubletabs_zinit = dblarr(ncols,nrows,comps) + zgal
    FOR I = 0, comps-1 DO doubletabs_zinit[*,*,I] += profileshifts[I]/1037.6167d
    doubletabs_siginit = dblarr(ncols,nrows,comps)
    FOR I = 0, N_ELEMENTS(profilesig)-1 DO $
        doubletabs_siginit[*,*,I] += profilesig[I]
    doubletabs_siglim = [1d,3000d]
    doubletabs_fix = bytarr(ncols,nrows,comps,4)
    doubletabs_cfinit = coveringfactor
    doubletabs_tauinit = opticaldepth

    init = {$
;
      plotindex:lineplotind,$
      fitindex:linefitind,$
      fcnfitdoublet: 'ifsf_doubletfcn',$
      fcninitpar: 'ifsf_initdoublet',$
;
      maxncomp: comps,$
;
      ndoubletabs: ndoubletabs,$
      doubletabs_cfinit: doubletabs_cfinit,$
      doubletabs_tauinit: doubletabs_tauinit,$
      doubletabs_zinit: doubletabs_zinit,$
      doubletabs_siginit: doubletabs_siginit,$
      doubletabs_siglim: doubletabs_siglim,$
      doubletabs_fix: doubletabs_fix,$
;
      ndoubletem: 0,$
;
      galaxy: gal, $
      zsys_gas: zgal,$
      outdir: directoryname+'/'+gal+'/',$
      wavelength: wavelength, $
      relativeflux: relativeflux, $
      error: relativeerror, $
      continuum: continuum, $
      flux: flux $
    }

    return,init

END
