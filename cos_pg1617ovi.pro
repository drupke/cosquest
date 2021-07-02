FUNCTION cos_pg1617ovi, directoryname, gal, zgal, profileshifts, $
                        profilesig, coveringfactor,opticaldepth

  fittedline = 'OVI'
  bad = 1d99
  gal = gal[0]
  ncols = 1
  nrows = 1
  centcol = 1
  centrow = 1
  zgal=zgal[0]
  comps=N_ELEMENTS(profileshifts)
  readcol,directoryname+'/'+gal+'/'+gal+'.txt', wavelength, flux, error, $
          FORMAT='D,D,D',/silent

  ; Finding the index to fit over
;  linefitreg=[1136,1152.2]
;  lineplotreg=[1132.01,1156]
;  contplotreg=[1132.01,1156]
  linefitreg=[1133.03d,1152.2d]
  lineplotreg=[1133.03d,1156d]
  contplotreg=[1133.03d,1175d]
  contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
     VALUE_LOCATE(wavelength,contplotreg[1])]
  linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
     VALUE_LOCATE(wavelength,linefitreg[1])]
  lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
     VALUE_LOCATE(wavelength,lineplotreg[1])]
  goodind = [[1141,1142],[1147.5,1151.5],[1153,1189],[1191,1192.5],[1194,1198],$
   [1201.5,1205],[1207.5,1224.5],[1225.5,1231.5],[1232.5,1236.5],$
   [1237.7,1238.2],[1239.2,1239.3],[1240.5,1242],[1243.5,1250],$
   [1251.5,1253.5],[1254.5,1257.5],[1261,1276.5],[1278.5,1301],$
   [1305,1320],[1321,1325],$
   [1326,1327],[1330,1333.5],[1336.5,1339],$
   [1344,1346.5],[1347.5,1349],[1350.5,1360],[1369,1369.8],[1374,1392.5],$
   [1394.5,1402],$
   [1404,1440]]
  for i=0,n_elements(goodind[0,*])-1 do begin
     newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
        VALUE_LOCATE(wavelength,goodind[0,i]),$
        START=VALUE_LOCATE(wavelength,goodind[0,i]))
     if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
  endfor
  weight=1d/error^2
  contfitreg=[[1133.03d,1175d]]
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
  twave *= 1d + zgal
  itflux = double(ifsf_interptemp(wavelength,twave,tflux))
  itflux /= median(itflux)

;  contfitreg=[[1132.01,1148],[1148,1156]]
;  fitfcn=['ifsf_fitpoly','ifsf_fitpoly']
;  fitargs=HASH()
;  fitargs['reg1'] = {fitord:3}
;  fitargs['reg2'] = {fitord:1}

  set_plot,'z'
  cgplot, wavelength, flux, XRAN=contplotreg, $
     YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
     1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
     XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
     xtit='Wavelength ($\Angstrom$)',$
     ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
;     continuum=ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
     continuum=ifsf_fitmulticont(wavelength, flux, weight, wavelength,itflux, $
     indextoplot,0,fitreg=contfitreg,$
     fitfcn=fitfcn, fitargs=fitargs)
  cgoplot, wavelength, continuum, color='Red',thick=4
  img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
     '_continuum',/jpeg,/nodialog,quality=100)

  relativeflux=flux/continuum
  relativeerror=error/continuum
  
; Removing problem areas by setting error to bad so that MPFITFUN ignores it
; intrinsic Lybeta; NI lines; FeII lines
;  relativeerror[VALUE_LOCATE(wavelength,1138):$
;                VALUE_LOCATE(wavelength,1139)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1133.8):$
                VALUE_LOCATE(wavelength,1135.1)] = bad
;  relativeerror[VALUE_LOCATE(wavelength,1142.7):$
;                VALUE_LOCATE(wavelength,1143.5)] = bad
;  relativeerror[VALUE_LOCATE(wavelength,1144.5):$
;                VALUE_LOCATE(wavelength,1145.2)] = bad


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Parameters for doublet fit
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ndoubletabs = dblarr(ncols,nrows)+comps
    doubletabs_zinit = dblarr(ncols,nrows,comps) + zgal
    FOR I = 0, comps-1 DO doubletabs_zinit[*,*,I] += profileshifts[I]/1037.6167d
    doubletabs_siginit = dblarr(ncols,nrows,comps)
    FOR I = 0, N_ELEMENTS(profilesig)-1 DO $
        doubletabs_siginit[*,*,I] += profilesig[I]
    doubletabs_siglim = [1d,1000d]
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
