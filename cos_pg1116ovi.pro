FUNCTION cos_pg1116ovi, directoryname, gal, zgal, profileshifts, $
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
  comps=N_ELEMENTS(profileshifts)
  readcol,directoryname+'/'+gal+'/'+gal+'.txt', wavelength, flux, error, $
          /silent,FORMAT='D,D,D'

  ; Finding the index to fit over
   linefitreg=[1195,1225]
   lineplotreg=[1195,1225]
   contplotreg=[1195,1225]
   contplotind=[VALUE_LOCATE(wavelength,contplotreg[0]),$
                VALUE_LOCATE(wavelength,contplotreg[1])]
   linefitind=[VALUE_LOCATE(wavelength,linefitreg[0]),$
               VALUE_LOCATE(wavelength,linefitreg[1])]
   lineplotind=[VALUE_LOCATE(wavelength,lineplotreg[0]),$
                VALUE_LOCATE(wavelength,lineplotreg[1])]
   goodind = [[1195,1195.5],[1197.5,1198.7],[1201.5,1202],[1203,1203.2],$
;              [1204.3,1205.7],[1207.8,1208.8],[1210.2,1210.5],[1211.2,1221],$
;              [1222.5,1225]]
              [1204.3,1205.7],[1207.8,1208.8],[1210.2,1210.5],[1211.2,1214],$
              [1217.5,1221],[1222.5,1225]]
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
   contfitreg=[[1195,1214],[1217.5,1225]]
   fitfcn=['ifsf_fitspline','ifsf_fitspline']
   fitargs=HASH()
   fitargs['reg1'] = {argsbkpts:{everyn:100}}
   fitargs['reg2'] = {argsbkpts:{everyn:100}}
  
   set_plot,'z'
   cgplot, wavelength, flux, XRAN=contplotreg, $
;      YRAN=[-.3*MAX(flux[contplotind[0]:contplotind[1]]),$
;      1.5*MAX(flux[contplotind[0]:contplotind[1]])],$
      YRAN=[0,1d-13],$
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

;  Removing problem areas by setting error to 0 so that MPFITFUN ignores it
   relativeerror[VALUE_LOCATE(wavelength,1199):$
                 VALUE_LOCATE(wavelength,1201.5)] = bad
   relativeerror[VALUE_LOCATE(wavelength,1206):$
                 VALUE_LOCATE(wavelength,1208)] = bad
   relativeerror[VALUE_LOCATE(wavelength,1214):$
                 VALUE_LOCATE(wavelength,1217.5)] = bad

   relativeflux[VALUE_LOCATE(wavelength,1214):$
                VALUE_LOCATE(wavelength,1217.5)] = 1d
                 
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
