FUNCTION cos_pg0953ovi, directoryname, gal, zgal, profileshifts, profilesig, $
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

; Finding the index to fit over
  linefitreg=[1268,1283]
  lineplotreg=[1263,1283]
  contplotreg=[1261,1283]
  contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                 VALUE_LOCATE(wavelength,contplotreg[1])]
  linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
              VALUE_LOCATE(wavelength,linefitreg[1])]
  lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
              VALUE_LOCATE(wavelength,lineplotreg[1])]
  goodind = [[1261,1264],[1265.5,1268.5],[1269.5,1270.4],[1271,1271.6],$
             [1272,1272.8],[1273.5,1275.7],[1276.5,1279.7],[1280.3,1283]]
  for i=0,n_elements(goodind[0,*])-1 do begin
     newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                   VALUE_LOCATE(wavelength,goodind[0,i]),$
                   START=VALUE_LOCATE(wavelength,goodind[0,i]))
     if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
  endfor
  weight=1d/error^2
  contfitreg=[[1261,1283]]
  fitfcn=strarr(n_elements(contfitreg[0,*]))+'ifsf_fitspline'
  fitargs=HASH()
  fitargs['reg1'] = {argsbkpts:{everyn:75}}

  set_plot,'z'
  cgplot, wavelength, flux, XRAN=contplotreg, $
          YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
          1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
          XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
          xtit='Wavelength ($\Angstrom$)',$
          ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
  continuum=ifsf_fitmulticont(wavelength, flux, weight, ignored,ignored, $
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
