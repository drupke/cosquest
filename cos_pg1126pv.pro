FUNCTION cos_pg1126pv, directoryname, gal, zgal, profileshifts, $
                       profilesig, coveringfactor,$
                       opticaldepth
                       
  fittedline = 'PV'
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
   linefitreg=[1175,1188]
   lineplotreg=[1170,1198]
   contplotreg=[1170,1198]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
               VALUE_LOCATE(wavelength,linefitreg[1])]
   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
                VALUE_LOCATE(wavelength,lineplotreg[1])]
   goodind = [[1170,1176],[1177,1186],[1187.8,1188.2],[1189.2,1189.8],$
              [1191.2,1192.7],[1194.5,1198]]
   for i=0,n_elements(goodind[0,*])-1 do begin
      newind=INDGEN(VALUE_LOCATE(wavelength,goodind[1,i])-$
                    VALUE_LOCATE(wavelength,goodind[0,i]),$
                    START=VALUE_LOCATE(wavelength,goodind[0,i]))
      if i eq 0 then indextoplot = newind else indextoplot = [indextoplot,newind]
   endfor
   weight = 1d/error^2
   contfitreg=[[1170,1198]]
   fitfcn=strarr(n_elements(contfitreg[0,*]))+'ifsf_fitspline'
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:140}}
    
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
  relativeerror[VALUE_LOCATE(wavelength,1190):$
                VALUE_LOCATE(wavelength,1191.5)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1192.5):$
                VALUE_LOCATE(wavelength,1194.5)] = bad
  relativeerror[VALUE_LOCATE(wavelength,1199):$
                VALUE_LOCATE(wavelength,1201.5)] = bad

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Parameters for doublet fit
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ndoubletabs = dblarr(ncols,nrows)+comps
    doubletabs_zinit = dblarr(ncols,nrows,comps) + zgal
    FOR I = 0, comps-1 DO doubletabs_zinit[*,*,I] += profileshifts[I]/1128.01d
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
