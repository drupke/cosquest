FUNCTION cos_pg1351nv, directoryname, gal, zgal, profileshifts, $
                       profilesig, coveringfactor,opticaldepth

  fittedline = 'NV'
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

; spline fit
; Finding the index to fit over
;  linefitreg=[1338,1350]
;  lineplotreg=[1336,1355]
;  contplotreg=[1330,1355]
;  contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
;     VALUE_LOCATE(wavelength,contplotreg[1])]
;  linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
;     VALUE_LOCATE(wavelength,linefitreg[1])]
;  lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
;     VALUE_LOCATE(wavelength,lineplotreg[1])]
;  goodind = [[1336,1338.5],[1341.5,1342],[1345.7,1346.5],[1350,1355]]
;  for i=0,n_elements(goodind[0,*])-1 do begin
;     newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
;        VALUE_LOCATE(wavelength,goodind[0,i]),$
;        START=VALUE_LOCATE(wavelength,goodind[0,i]))
;     if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
;  endfor
;  weight=1d/error^2
;  contfitreg=[[1336,1355]]
;  fitfcn=strarr(n_elements(contfitreg[0,*]))+'ifsf_fitspline'
;  fitargs=HASH()
;  fitargs['reg1'] = {argsbkpts:{everyn:90}}
  
; template fit
  linefitreg=[1338,1350]
  lineplotreg=[1330,1355]
  contplotreg=[1325,1365]
  contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
     VALUE_LOCATE(wavelength,contplotreg[1])]
  linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
     VALUE_LOCATE(wavelength,linefitreg[1])]
  lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
     VALUE_LOCATE(wavelength,lineplotreg[1])]
  goodind = [[1330,1333],[1336,1338.5],[1345.7,1346.5],[1350,1360]]
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
   contfitreg=[[1330,1360]]
   fitfcn=['ifsf_fittemplate']
   fitargs=HASH()
   parinfo = replicate({value:0d,fixed:0b},4)
   parinfo[0].value = 1d-14
;   parinfo[3].fixed = 1b
   fitargs['reg1'] = {templatefcn: 'cos_quasarcomposite',$
                      parinfo:parinfo, npar: 4}
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
          YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
                1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
          XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
          xtit='Wavelength ($\Angstrom$)',$
          ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
;  continuum=ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
  continuum=ifsf_fitmulticont(wavelength, flux, weight, wavelength, itflux, $
                              indextoplot,0,fitreg=contfitreg,$
                              fitfcn=fitfcn, fitargs=fitargs)
  cgoplot, wavelength, continuum, color='Red',thick=4
  img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
                   '_continuum',/jpeg,/nodialog,quality=100)
  
  relativeflux=flux/continuum
  relativeerror=error/continuum


  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Parameters for doublet fit
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ndoubletabs = dblarr(ncols,nrows)+comps
    doubletabs_zinit = dblarr(ncols,nrows,comps) + zgal
    FOR I = 0, comps-1 DO doubletabs_zinit[*,*,I] += profileshifts[I]/1242.804d
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
      argslinelist: {vacuum:1b},$
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
         cfcorr: 0b,$
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
