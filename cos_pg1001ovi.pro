FUNCTION cos_pg1001ovi, directoryname, gal, zgal, profileshifts, profilesig, $
                        coveringfactor,opticaldepth

  fittedline = 'OVI'
  bad = 1d99
  gal = gal[0]
  ncols = 1
  nrows = 1
  zgal=zgal[0]
  comps=N_ELEMENTS(profileshifts)
  readcol,directoryname+'/'+gal+'/'+gal+'.txt', wavelength, flux, error, $
          FORMAT='D,D,D',/silent

; poly/bspline fit
;  Finding the index to fit over
  linefitreg=[1174,1206]
  lineplotreg=[1160,1210]
  contplotreg=[1130,1210]
  contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                 VALUE_LOCATE(wavelength,contplotreg[1])]
  linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
              VALUE_LOCATE(wavelength,linefitreg[1])]
  lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
              VALUE_LOCATE(wavelength,lineplotreg[1])]
goodind = [[1130,1133],[1136,1144],[1146,1152],[1153.5,1160],$
   [1188,1190],[1191,1193],$
   [1194.5,1197.5],[1199,1199.2],[1199.8,1200],$
   [1201.2,1203.5],$
   [1205.7,1206],[1207.7,1210]]
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
;  contfitreg=[[1194.5,1210]]
  contfitreg=[[1130,1157],[1157,1192],[1192,1210]]
  fitfcn=['ifsf_fitspline','ifsf_fitpoly','ifsf_fitspline']
  fitargs=HASH()
  fitargs['reg1'] = {argsbkpts:{everyn:100}}
  fitargs['reg2'] = {}
  fitargs['reg3'] = {argsbkpts:{everyn:100}}

  set_plot,'z'
  cgplot, wavelength, flux, XRAN=contplotreg, $
     YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
     1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
     XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
     xtit='Wavelength ($\Angstrom$)',$
     ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
     continuum=ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
                                 indextoplot,0,fitreg=contfitreg,$
                                 fitfcn=fitfcn, fitargs=fitargs,quiet=1b)
  cgoplot, wavelength, continuum, color='Red',thick=4
  img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
     '_continuum',/jpeg,/nodialog,quality=100)

;;  template fit
;  linefitreg=[1174,1206]
;  lineplotreg=[1130,1210]
;  contplotreg=[1130,1210]
;;  contplotreg=[1194.5,1210]
;  contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;                 VALUE_LOCATE(wavelength,contplotreg[1])]
;  linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
;              VALUE_LOCATE(wavelength,linefitreg[1])]
;  lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
;              VALUE_LOCATE(wavelength,lineplotreg[1])]
;   goodind = [[1130,1133],[1136,1144],[1146,1152],[1153.5,1160],$
;             [1188,1190],[1191,1193],$
;             [1194.5,1197.5],[1199,1199.2],[1199.8,1200],$
;             [1201.2,1203.5],$
;             [1205.7,1206],[1207.7,1210]]
;  for i=0,n_elements(goodind[0,*])-1 do begin
;     newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
;                   VALUE_LOCATE(wavelength,goodind[0,i]),$
;                   START=VALUE_LOCATE(wavelength,goodind[0,i]))
;     if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
;  endfor
;   weight=dblarr(n_elements(error))
;   tolerance = 1d-80
;   ibd = where(error^2d LT tolerance)
;   igd = where(error^2d GE tolerance)
;   weight[ibd] = 0.0d0
;   weight[igd] = 1d/error[igd]^2d
;   contfitreg=[[1130,1210]]
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
;   twave *= 1d + 0.16d
;   itflux = double(ifsf_interptemp(wavelength,twave,tflux))
;   itflux /= median(itflux)
;   
;   set_plot,'z'
;   cgplot, wavelength, flux, XRAN=contplotreg, $
;;           YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;;                 1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
;           YRAN=[0,3d-14],$
;           XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
;           xtit='Wavelength ($\Angstrom$)',$
;           ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
;   continuum=$
;      ifsf_fitmulticont(wavelength, flux, weight, wavelength, itflux, $
;                        indextoplot,0,fitreg=contfitreg,$
;                        fitfcn=fitfcn, fitargs=fitargs,quiet=1b)
;   cgoplot, wavelength, continuum, color='Red',thick=4
;   img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
;                    '_continuum',/jpeg,/nodialog,quality=100)


  relativeflux=flux/continuum
  relativeerror=error/continuum

  ; Removing problem areas by setting error to 0 so that MPFITFUN ignores it
  relativeerror[VALUE_LOCATE(wavelength,1189.5):$
    VALUE_LOCATE(wavelength,1191.5)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1192.5):$
    VALUE_LOCATE(wavelength,1194.5)] = bad

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Parameters for doublet fit
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ndoubletabs = dblarr(ncols,nrows)+comps
    doubletabs_zinit = dblarr(ncols,nrows,comps) + zgal
    FOR I = 0, comps-1 DO doubletabs_zinit[*,*,I] += profileshifts[I]/1037.6167d
    doubletabs_siginit = dblarr(ncols,nrows,comps)
    FOR I = 0, N_ELEMENTS(profilesig)-1 DO $
        doubletabs_siginit[*,*,I] += profilesig[I]
    doubletabs_siglim = [1d,2000d]
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
