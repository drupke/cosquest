FUNCTION cos_pg1126nv, directoryname, gal, zgal, profileshifts, $
                       profilesig, coveringfactor,$
                       opticaldepth
                       
  fittedline = 'NV'
  bad = 1d99
  gal = gal[0]
  bin = 2d
  ncols = 1
  nrows = 1
  centcol = 1
  centrow = 1
  zgal=zgal[0]
  comps=N_ELEMENTS(profileshifts)
  readcol,directoryname+'/'+gal+'/'+gal+'.txt', wavelength, flux, error, $
          FORMAT='D,D,D',/silent

  ; Finding the index to fit over
   linefitreg=[1292,1319]
   lineplotreg=[1292,1319]
;   lineplotreg=[1295,1308]
   contplotreg=[1290,1328]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
               VALUE_LOCATE(wavelength,linefitreg[1])]
   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
                VALUE_LOCATE(wavelength,lineplotreg[1])]
   goodind = [[1290,1292],[1296.5,1297.5],[1301.5,1301.7],$
              [1306.3,1306.7],[1310.6,1310.95],[1314,1315.2],[1315.7,1316.8],$
              [1318,1328]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
   endfor
   weight = 1d/error^2
   contfitreg=[[1290,1314.5],[1314.5,1328]]
   fitfcn=strarr(n_elements(contfitreg[0,*]))+'ifsf_fitspline'
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:140}}
   fitargs['reg2'] = {argsbkpts:{everyn:100}}
;   contfitreg=[[1290,1297],[1297,1306.5],[1306.5,1314.5],[1314.5,1328]]
;   fitfcn=['ifsf_fitspline','ifsf_fitpoly','ifsf_fitpoly','ifsf_fitspline']
;   fitargs=HASH()
;   fitargs['reg1'] = {argsbkpts:{everyn:100}}
;   fitargs['reg2'] = {fitord:1}
;   fitargs['reg3'] = {fitord:1}
;   fitargs['reg4'] = {argsbkpts:{everyn:100}}
    
  set_plot,'z'
  cgplot, wavelength, flux, XRAN=contplotreg, $
          YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
                1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
          XSTYLE=1,YSTYLE=1,backg='Black',axiscolor='White',color='White',$
          xtit='Wavelength ($\Angstrom$)',$
          ytit='Flux (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)'
  continuum=ifsf_fitmulticont(wavelength, flux, weight, ignored, ignored, $
                              indextoplot,0,fitreg=contfitreg,$
                              fitfcn=fitfcn, fitargs=fitargs)
  cgoplot, wavelength, continuum, color='Red',thick=4
  img = cgsnapshot(filename=directoryname+'/'+gal+'/'+'/'+gal+fittedline+$
                   '_continuum',/jpeg,/nodialog,quality=100)
  
  relativeflux=flux/continuum
  relativeerror=error/continuum
  
  ; Removing problem areas by setting error to something big
  ; so that MPFITFUN ignores it
  relativeerror[VALUE_LOCATE(wavelength, 1301.8):$
                VALUE_LOCATE(wavelength,1302.6)] = 1d99

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
