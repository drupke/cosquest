; docformat = 'rst'
;
;+
;
; Plot correlations.
;
; :Categories:
;    COSQUEST
;
; :Returns:
;    Postscript plots.
;
; :Params:
;
; :Keywords:
;
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2016aug02, DSNR, created
;
; :Copyright:
;    Copyright (C) 2016--2021 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
; Routine from E. Cameron 2011, PASP
; z = DINDGEN(10000)*0.0001d
function betaprob,c,k,n,z

   Beta = IBETA(k+1,n-k+1,z)
   ill = VALUE_LOCATE(Beta,(1-c)/2)
   iul = VALUE_LOCATE(Beta,1-(1-c)/2)
   plower = z[ill]
   pupper = z[iul]

   return,[plower,pupper]

end
;;;
pro cos_corr

   bad = 1d99
   fwhm2sig = 2d*sqrt(2d*alog(2d))
   lsun_ergps = 3.839d33 ; L_Sun in erg/s
   llsun_ergps = alog10(lsun_ergps)
   c_cms = 2.99792458d10
   c_kms = 2.99792458d5
   spy = 3600d*24d*365.25d
   msun_g = 1.989d33
   eddfac = 1.26d38 ; 4 Pi G m_p c / sigma_T, in erg/s/M_sun

   ; linmix_err iterations
   lm_miniter = 7894L  ; 15787 is 4sigma, and miniter typically gets doubled
   lm_maxiter = 315740L  ; 20 x 15787

   ; Threshold pval for marking as significant in output regression table
   tpval = 0.05

   qsotab='/Users/drupke/Box Sync/qsos/qsos.csv'
   cosdir='/Users/drupke/Box Sync/cosquest/'
   fitdir=cosdir+'fits/'
   plotdir=cosdir+'plots/correlations/'
   tabdir=cosdir+'tables/'
   specdir=cosdir+'spectra/'

   linelist = ifsf_linelist(['OVI1031','OVI1037','NV1238','NV1242',$
      'PV1117','PV1128'],vacuum=1b)

;  for tables
   lerr0 = '$_{-'
   lerr1 = '}^{+'
   lerr2 = '}$'
   amp = '&'
   dslash = '\\'
   ndat = '\nodata'
   lineofdashes = strjoin(replicate('-',62))


;  Read table
   rows=[3,85]
   gal = read_csvcol(qsotab,'A',rows=rows,sep=',',type='string')
   sgal = read_csvcol(qsotab,'C',rows=rows,sep=',',type='string')
   z = read_csvcol(qsotab,'D',rows=rows,sep=',',junk=bad)
   lir = read_csvcol(qsotab,'F',rows=rows,sep=',',junk=bad)
   lbol = read_csvcol(qsotab,'G',rows=rows,sep=',',junk=bad)
   agnfrac = read_csvcol(qsotab,'H',rows=rows,sep=',',junk=bad)
   agnfraclb = read_csvcol(qsotab,'I',rows=rows,sep=',',junk=bad)
   agnfracub = read_csvcol(qsotab,'J',rows=rows,sep=',',junk=bad)
   mbh_phot = read_csvcol(qsotab,'M',rows=rows,sep=',',junk=bad)
   lmbh_rev = read_csvcol(qsotab,'N',rows=rows,sep=',',junk=bad)
   lmbh_rev_errhi = read_csvcol(qsotab,'O',rows=rows,sep=',',junk=bad)
   lmbh_rev_errlo = read_csvcol(qsotab,'P',rows=rows,sep=',',junk=bad)
   lmbh_rev_type = read_csvcol(qsotab,'Q',rows=rows,sep=',',type='string')
   cossamp = read_csvcol(qsotab,'AD',rows=rows,sep=',',type='string')
   nvstatus = read_csvcol(qsotab,'AE',rows=rows,sep=',',type='string')
   ovistatus = read_csvcol(qsotab,'AF',rows=rows,sep=',',type='string')
   gamxray = read_csvcol(qsotab,'AR',rows=rows,sep=',',junk=bad)
   gamxray_errlo = read_csvcol(qsotab,'AS',rows=rows,sep=',',junk=bad)
   gamxray_errhi = read_csvcol(qsotab,'AT',rows=rows,sep=',',junk=bad)
   nhxray = read_csvcol(qsotab,'AU',rows=rows,sep=',',junk=bad)
   nhxray_errlo = read_csvcol(qsotab,'AV',rows=rows,sep=',',junk=bad)
   nhxray_errhi = read_csvcol(qsotab,'AW',rows=rows,sep=',',junk=bad)
   fsoftxray = read_csvcol(qsotab,'AX',rows=rows,sep=',',junk=bad)
   fsoftxray_errlo = read_csvcol(qsotab,'AY',rows=rows,sep=',',junk=bad)
   fsoftxray_errhi = read_csvcol(qsotab,'AZ',rows=rows,sep=',',junk=bad)
   fhardxray = read_csvcol(qsotab,'BA',rows=rows,sep=',',junk=bad)
   fhardxray_errlo = read_csvcol(qsotab,'BB',rows=rows,sep=',',junk=bad)
   fhardxray_errhi = read_csvcol(qsotab,'BC',rows=rows,sep=',',junk=bad)
   lsoftxray = read_csvcol(qsotab,'BD',rows=rows,sep=',',junk=bad)
   lhardxray = read_csvcol(qsotab,'BE',rows=rows,sep=',',junk=bad)
   alphaox = read_csvcol(qsotab,'BL',rows=rows,sep=',',junk=bad)
   lfir = read_csvcol(qsotab,'BN',rows=rows,sep=',',junk=bad)
   ebvgal = read_csvcol(qsotab,'BS',rows=rows,sep=',',junk=bad)

;  Get flam1125
   readcol,tabdir+'tab_specvals.txt',spgal,flam1125,skip=4,/silent,format='(A,D,X)'

;  Parse table data
   icos = where(cossamp eq 'V18',ncos)
   nlines = n_elements(gal)
   tabdat = orderedhash()
   tabdat['gal'] = gal[icos]
   tabdat['sgal'] = sgal[icos]
   tabdat['z'] = z[icos]
   tabdat['lir'] = lir[icos]
   tabdat['lfir'] = lfir[icos]
   tabdat['lbol'] = lbol[icos]
   tabdat['agnfrac'] = agnfrac[icos]
   tabdat['agnfraclb'] = agnfraclb[icos]
   tabdat['agnfracub'] = agnfracub[icos]
   tabdat['mbh_phot'] = mbh_phot[icos]
   tabdat['lmbh_rev'] = lmbh_rev[icos]
   tabdat['lmbh_rev_errlo'] = lmbh_rev_errlo[icos]
   tabdat['lmbh_rev_errhi'] = lmbh_rev_errhi[icos]
   tabdat['lmbh_rev_type'] = lmbh_rev_type[icos]
   tabdat['nvstatus'] = nvstatus[icos]
   tabdat['ovistatus'] = ovistatus[icos]
   tabdat['alphaox'] = alphaox[icos]
   tabdat['ebvgal'] = ebvgal[icos]
   tabdat['flam1125'] = flam1125

;  Correct cosmologies
;  To convert cosmologies, back out -2 log(d_L)_old and add 2 log(d_L)_new
   tabdat['lbol'] = tabdat['lbol'] - $
      2d*(alog10(lumdist(tabdat['z'],H0=70d,Om=0.3d,Lam=0.7d)) - $
      alog10(lumdist(tabdat['z'],H0=69.3d,Om=0.287d,Lam=0.713d)))
   tabdat['lir'] = tabdat['lir'] - $
      2d*(alog10(lumdist(tabdat['z'],H0=67.8d,Om=0.308d,Lam=0.692d)) - $
      alog10(lumdist(tabdat['z'],H0=69.3d,Om=0.287d,Lam=0.713d)))

;  output table for x-ray quantities for Table 1
   openw,lun_tmp,tabdir+'tab_table1quants_xray.tex',/get_lun
   printf,lun_tmp,'#Col 2: log[L(0.5-2keV)/10^44 erg/s]_error_low^error_hi
   printf,lun_tmp,'#Col 3: log[L(2-10keV)/10^44 erg/s]_error_low^error_hi
   
;  X-ray data requires special handling because of multiple measurements
   tabdat['n_xray'] = dblarr(ncos)
   tabdat['gamxray'] = dblarr(ncos,6)+bad
   tabdat['gamxray_errlo'] = dblarr(ncos,6)
   tabdat['gamxray_errhi'] = dblarr(ncos,6)
   tabdat['nhxray'] = dblarr(ncos,6)+bad
   tabdat['nhxray_errlo'] = dblarr(ncos,6)
   tabdat['nhxray_errhi'] = dblarr(ncos,6)
   tabdat['nhxray_lim'] = dblarr(ncos,6)+bad
   tabdat['fsoftxray'] = dblarr(ncos,6)+bad
   tabdat['fsoftxray_errlo'] = dblarr(ncos,6)
   tabdat['fsoftxray_errhi'] = dblarr(ncos,6)
   tabdat['fhardxray'] = dblarr(ncos,6)+bad
   tabdat['fhardxray_errlo'] = dblarr(ncos,6)
   tabdat['fhardxray_errhi'] = dblarr(ncos,6)
   tabdat['lsoftxray'] = dblarr(ncos,6)+bad
   tabdat['lsoftxray_errlo'] = dblarr(ncos,6)
   tabdat['lsoftxray_errhi'] = dblarr(ncos,6)
   tabdat['lhardxray'] = dblarr(ncos,6)+bad
   tabdat['lhardxray_errlo'] = dblarr(ncos,6)
   tabdat['lhardxray_errhi'] = dblarr(ncos,6)
   tabdat['ltotxray'] = dblarr(ncos,6)+bad
   tabdat['ltotxray_errlo'] = dblarr(ncos,6)
   tabdat['ltotxray_errhi'] = dblarr(ncos,6)
   tabdat['lsoftratxray'] = dblarr(ncos,6)+bad
   tabdat['lsoftratxray_err'] = dblarr(ncos,6)+bad
   tabdat['lxlbol'] = dblarr(ncos,6)+bad
   tabdat['lxlbol_errlo'] = dblarr(ncos,6)+bad
   tabdat['lxlbol_errhi'] = dblarr(ncos,6)+bad
   j = -1 ; galaxy index
   k = 0 ; zero x-ray component index
   for i=0,nlines-1 do begin
      if gal[i] ne '' AND cossamp[i] eq 'V18' then begin
         j++ ; increment galaxy index
         hardstring=''
         softstring=''
         if gamxray[i] ne bad then begin
            k = 0 ; re-zero x-ray component index
;           record first x-ray measurement
            tabdat['n_xray',j] = k+1
            tabdat['gamxray',j,k] = gamxray[i]
            tabdat['gamxray_errlo',j,k] = gamxray_errlo[i]
            tabdat['gamxray_errhi',j,k] = gamxray_errhi[i]
            if nhxray[i] eq 0d then begin
                if nhxray_errhi[i] ne 0d then $
                   tabdat['nhxray_lim',j,k] = alog10(nhxray_errhi[i]) + 22d $
                else $
                   tabdat['nhxray_lim',j,k] = 21d
            endif else begin
               tabdat['nhxray',j,k] = alog10(nhxray[i])+22d
               tabdat['nhxray_errlo',j,k] = $
                  alog10(nhxray[i]) - alog10(nhxray[i] - nhxray_errlo[i])
               tabdat['nhxray_errhi',j,k] = $
                  alog10(nhxray[i] + nhxray_errhi[i]) - alog10(nhxray[i])
            endelse
            if lsoftxray[i] ne bad then begin
               tabdat['lsoftxray',j,k] = alog10(lsoftxray[i])+44d
               softstring = $
                  string(tabdat['lsoftxray',j,k],format='(D0.2)')
               if fsoftxray[i] ne bad then begin
                  tabdat['fsoftxray',j,k] = alog10(fsoftxray[i]) - 12d
                  if fsoftxray_errlo[i] ne 0d AND $
                     fsoftxray_errlo[i] ne bad then begin
                     ; compute flux errors in log space
                     tabdat['fsoftxray_errlo',j,k] = $
                        alog10(fsoftxray[i]) - $
                        alog10(fsoftxray[i] - fsoftxray_errlo[i])
                     tabdat['fsoftxray_errhi',j,k] = $
                        alog10(fsoftxray[i] + fsoftxray_errhi[i]) - $
                        alog10(fsoftxray[i])
                     fsoftxray_err = (fsoftxray_errlo[i] + fsoftxray_errhi[i])/2d
                     ; Compute luminosity errors from flux errors
                     tabdat['lsoftxray_errlo',j,k] = $
                        alog10(lsoftxray[i]) - $
                        alog10(lsoftxray[i] - lsoftxray[i] * $
                        fsoftxray_errlo[i] / fsoftxray[i])
                     tabdat['lsoftxray_errhi',j,k] = $
                        alog10(lsoftxray[i] + lsoftxray[i] * $
                        fsoftxray_errhi[i] / fsoftxray[i]) - $
                        alog10(lsoftxray[i])
                     if tabdat['lsoftxray_errlo',j,k] lt 0.0095 OR $
                        tabdat['lsoftxray_errhi',j,k] lt 0.0095 then $
                        fstring = '(D0.3,A0,D0.3,A0,D0.3,A0)' $
                     else fstring = '(D0.2,A0,D0.2,A0,D0.2,A0)'
                     softstring = $
                        string(tabdat['lsoftxray',j,k],'$_{-',$
                        tabdat['lsoftxray_errlo',j,k],'}^{+',$
                        tabdat['lsoftxray_errhi',j,k],'}$',$
                        format=fstring)
                   ; case of no flux errors; assume 10%
                  endif else begin
                     tabdat['fsoftxray_errlo',j,k] = $
                        alog10(fsoftxray[i]) -  alog10(0.9d * fsoftxray[i])
                     tabdat['fsoftxray_errhi',j,k] = $
                        alog10(1.1d * fsoftxray[i]) - alog10(fsoftxray[i])
                     tabdat['lsoftxray_errlo',j,k] = $
                        alog10(lsoftxray[i]) -  alog10(0.9d * lsoftxray[i])
                     tabdat['lsoftxray_errhi',j,k] = $
                        alog10(1.1d * lsoftxray[i]) - alog10(lsoftxray[i])
                     fsoftxray_err = 0.1d * fsoftxray[i]
                  endelse
               endif
            endif
            if lhardxray[i] ne bad then begin
               tabdat['lhardxray',j,k] = alog10(lhardxray[i])+44d
               hardstring = $
                  string(tabdat['lhardxray',j,k],format='(D0.2)')
               if fhardxray[i] ne bad then begin
                  tabdat['fhardxray',j,k] = alog10(fhardxray[i]) - 12d
                  if fhardxray_errlo[i] ne 0d AND $
                     fhardxray_errlo[i] ne bad then begin
                     ; compute flux errors in log space
                     tabdat['fhardxray_errlo',j,k] = $
                        alog10(fhardxray[i]) - $
                        alog10(fhardxray[i] - fhardxray_errlo[i])
                     tabdat['fhardxray_errhi',j,k] = $
                        alog10(fhardxray[i] + fhardxray_errhi[i]) - $
                        alog10(fhardxray[i])
                     tabdat['lhardxray_errlo',j,k] = $
                        alog10(lhardxray[i]) - $
                        alog10(lhardxray[i] - lhardxray[i] * $
                        fhardxray_errlo[i] / fhardxray[i])
                     fhardxray_err = (fhardxray_errlo[i] + fhardxray_errhi[i])/2d
                     ; Compute luminosity errors from flux errors
                     tabdat['lhardxray_errhi',j,k] = $
                        alog10(lhardxray[i] + lhardxray[i] * $
                        fhardxray_errhi[i] / fhardxray[i]) - $
                        alog10(lhardxray[i])
                     if tabdat['lhardxray_errlo',j,k] lt 0.0095 OR $
                        tabdat['lhardxray_errhi',j,k] lt 0.0095 then $
                        fstring = '(D0.3,A0,D0.3,A0,D0.3,A0)' $
                     else fstring = '(D0.2,A0,D0.2,A0,D0.2,A0)'
                     hardstring = $
                        string(tabdat['lhardxray',j,k],'$_{-',$
                        tabdat['lhardxray_errlo',j,k],'}^{+',$
                        tabdat['lhardxray_errhi',j,k],'}$',$
                        format=fstring)
                   ; case of no flux errors; assume 10%
                  endif else begin
                     tabdat['fhardxray_errlo',j,k] = $
                        alog10(fhardxray[i]) -  alog10(0.9d * fhardxray[i])
                     tabdat['fhardxray_errhi',j,k] = $
                        alog10(1.1d * fhardxray[i]) - alog10(fhardxray[i])
                     tabdat['lhardxray_errlo',j,k] = $
                        alog10(lhardxray[i]) -  alog10(0.9d * lhardxray[i])
                     tabdat['lhardxray_errhi',j,k] = $
                        alog10(1.1d * lhardxray[i]) - alog10(lhardxray[i])
                     fhardxray_err = 0.1d * fhardxray[i]
                  endelse
               endif
            endif
            if lhardxray[i] ne bad AND lsoftxray[i] ne bad AND $
               lsoftxray[i] ne 0d then begin
               tabdat['ltotxray',j,k] = alog10(lsoftxray[i] + lhardxray[i])+44d
               totxray = lsoftxray[i] + lhardxray[i]
               tabdat['ltotxray_errlo',j,k] = $
                  sqrt(tabdat['lsoftxray_errlo',j,k]^2d + $
                     tabdat['lhardxray_errlo',j,k]^2d)
               tabdat['ltotxray_errhi',j,k] = $
                  sqrt(tabdat['lsoftxray_errhi',j,k]^2d + $
                     tabdat['lhardxray_errhi',j,k]^2d)
               tabdat['lsoftratxray',j,k] = $
                  tabdat['lsoftxray',j,k] - tabdat['ltotxray',j,k]
               tabdat['lsoftratxray_err',j,k] = $
                     sqrt(fsoftxray[i]^2d * fhardxray_err^2d + $
                     fhardxray[i]^2d * fsoftxray_err^2d) / $
                     (fsoftxray[i]+fhardxray[i])^2d
               tabdat['lxlbol',j,k] = tabdat['ltotxray',j,k] - $
                  llsun_ergps - tabdat['lbol',j]
               tabdat['lxlbol_errlo',j,k] = tabdat['ltotxray_errlo',j,k]
               tabdat['lxlbol_errhi',j,k] = tabdat['ltotxray_errhi',j,k]
            endif
            printf,lun_tmp,gal[i],softstring,hardstring,$
                format='(A-12,2A27)'
            k++
         endif else begin
            tabdat['n_xray',j] = 0 ; record number of x-ray components
            k = 0 ; re-zero x-ray component index
         endelse
      endif else if gal[i] eq '' then begin
         if k gt 0 then begin
            hardstring=''
            softstring=''
            if gamxray[i] ne bad then begin
               tabdat['n_xray',j] = k+1
               tabdat['gamxray',j,k] = gamxray[i]
               tabdat['gamxray_errlo',j,k] = gamxray_errlo[i]
               tabdat['gamxray_errhi',j,k] = gamxray_errhi[i]
               if nhxray[i] eq 0d then begin
                  if nhxray_errhi[i] ne 0d then $
                     tabdat['nhxray_lim',j,k] = alog10(nhxray_errhi[i]) + 22d $
                  else $
                     tabdat['nhxray_lim',j,k] = 21d
               endif else begin
                  tabdat['nhxray',j,k] = alog10(nhxray[i])+22d
                  tabdat['nhxray_errlo',j,k] = $
                     alog10(nhxray[i]) - alog10(nhxray[i] - nhxray_errlo[i])
                  tabdat['nhxray_errhi',j,k] = $
                     alog10(nhxray[i] + nhxray_errhi[i]) - alog10(nhxray[i])
               endelse
               if lsoftxray[i] ne bad then begin
                  tabdat['lsoftxray',j,k] = alog10(lsoftxray[i])+44d
                  softstring = $
                     string(tabdat['lsoftxray',j,k],format='(D0.2)')
                  if fsoftxray[i] ne bad then begin
                     tabdat['fsoftxray',j,k] = alog10(fsoftxray[i]) - 12d
                     if fsoftxray_errlo[i] ne 0d AND $
                        fsoftxray_errlo[i] ne bad then begin
                        ; compute flux errors in log space
                        tabdat['fsoftxray_errlo',j,k] = $
                           alog10(fsoftxray[i]) - $
                           alog10(fsoftxray[i] - fsoftxray_errlo[i])
                        tabdat['fsoftxray_errhi',j,k] = $
                           alog10(fsoftxray[i] + fsoftxray_errhi[i]) - $
                           alog10(fsoftxray[i])
                        fsoftxray_err = (fsoftxray_errlo[i] + fsoftxray_errhi[i])/2d
                        ; Compute luminosity errors from flux errors
                        tabdat['lsoftxray_errlo',j,k] = $
                           alog10(lsoftxray[i]) - $
                           alog10(lsoftxray[i] - lsoftxray[i] * $
                           fsoftxray_errlo[i] / fsoftxray[i])
                        tabdat['lsoftxray_errhi',j,k] = $
                           alog10(lsoftxray[i] + lsoftxray[i] * $
                           fsoftxray_errhi[i] / fsoftxray[i]) - $
                           alog10(lsoftxray[i])
                        if tabdat['lsoftxray_errlo',j,k] lt 0.0095 OR $
                           tabdat['lsoftxray_errhi',j,k] lt 0.0095 then $
                           fstring = '(D0.3,A0,D0.3,A0,D0.3,A0)' $
                        else fstring ='(D0.2,A0,D0.2,A0,D0.2,A0)'
                        softstring = $
                           string(tabdat['lsoftxray',j,k],'$_{-',$
                           tabdat['lsoftxray_errlo',j,k],'}^{+',$
                           tabdat['lsoftxray_errhi',j,k],'}$',$
                           format=fstring)
                        ; case of no flux errors; assume 10%
                     endif else begin
                        tabdat['fsoftxray_errlo',j,k] = $
                           alog10(fsoftxray[i]) -  alog10(0.9d * fsoftxray[i])
                        tabdat['fsoftxray_errhi',j,k] = $
                           alog10(1.1d * fsoftxray[i]) - alog10(fsoftxray[i])
                        tabdat['lsoftxray_errlo',j,k] = $
                           alog10(lsoftxray[i]) -  alog10(0.9d * lsoftxray[i])
                        tabdat['lsoftxray_errhi',j,k] = $
                           alog10(1.1d * lsoftxray[i]) - alog10(lsoftxray[i])
                        fsoftxray_err = 0.1d * fsoftxray[i]
                     endelse
                  endif
               endif
               if lhardxray[i] ne bad then begin
                  tabdat['lhardxray',j,k] = alog10(lhardxray[i])+44d
                  hardstring = $
                     string(tabdat['lhardxray',j,k],format='(D0.2)')
                  if fhardxray[i] ne bad then begin
                     tabdat['fhardxray',j,k] = alog10(fhardxray[i]) - 12d
                     if fhardxray_errlo[i] ne 0d AND $
                        fhardxray_errlo[i] ne bad then begin
                        ; compute flux errors in log space
                        tabdat['fhardxray_errlo',j,k] = $
                           alog10(fhardxray[i]) - $
                           alog10(fhardxray[i] - fhardxray_errlo[i])
                        tabdat['fhardxray_errhi',j,k] = $
                           alog10(fhardxray[i] + fhardxray_errhi[i]) - $
                           alog10(fhardxray[i])
                        tabdat['lhardxray_errlo',j,k] = $
                           alog10(lhardxray[i]) - $
                           alog10(lhardxray[i] - lhardxray[i] * $
                           fhardxray_errlo[i] / fhardxray[i])
                        fhardxray_err = (fhardxray_errlo[i] + fhardxray_errhi[i])/2d
                        ; Compute luminosity errors from flux errors
                        tabdat['lhardxray_errhi',j,k] = $
                           alog10(lhardxray[i] + lhardxray[i] * $
                           fhardxray_errhi[i] / fhardxray[i]) - $
                           alog10(lhardxray[i])
                        if tabdat['lhardxray_errlo',j,k] lt 0.0095 OR $
                           tabdat['lhardxray_errhi',j,k] lt 0.0095 then $
                           fstring ='(D0.3,A0,D0.3,A0,D0.3,A0)' $
                        else fstring = '(D0.2,A0,D0.2,A0,D0.2,A0)'
                        hardstring = $
                           string(tabdat['lhardxray',j,k],'$_{-',$
                           tabdat['lhardxray_errlo',j,k],'}^{+',$
                           tabdat['lhardxray_errhi',j,k],'}$',$
                           format=fstring)
                        ; case of no flux errors; assume 10%
                     endif else begin
                        tabdat['fhardxray_errlo',j,k] = $
                           alog10(fhardxray[i]) -  alog10(0.9d * fhardxray[i])
                        tabdat['fhardxray_errhi',j,k] = $
                           alog10(1.1d * fhardxray[i]) - alog10(fhardxray[i])
                        tabdat['lhardxray_errlo',j,k] = $
                           alog10(lhardxray[i]) -  alog10(0.9d * lhardxray[i])
                        tabdat['lhardxray_errhi',j,k] = $
                           alog10(1.1d * lhardxray[i]) - alog10(lhardxray[i])
                        fhardxray_err = 0.1d * fhardxray[i]
                     endelse
                  endif
               endif
               if lhardxray[i] ne bad AND lsoftxray[i] ne bad $
                  AND lsoftxray[i] ne 0d then begin
                  tabdat['ltotxray',j,k] = alog10(lsoftxray[i] + lhardxray[i])+44d
                  totxray = lsoftxray[i] + lhardxray[i]
                  tabdat['ltotxray_errlo',j,k] = $
                     sqrt(tabdat['lsoftxray_errlo',j,k]^2d + $
                     tabdat['lhardxray_errlo',j,k]^2d)
                  tabdat['ltotxray_errhi',j,k] = $
                     sqrt(tabdat['lsoftxray_errhi',j,k]^2d + $
                     tabdat['lhardxray_errhi',j,k]^2d)
                  tabdat['lsoftratxray',j,k] = $
                     tabdat['lsoftxray',j,k] - tabdat['ltotxray',j,k]
                  tabdat['lsoftratxray_err',j,k] = $
                     sqrt(fsoftxray[i]^2d * fhardxray_err^2d + $
                     fhardxray[i]^2d * fsoftxray_err^2d) / $
                     (fsoftxray[i]+fhardxray[i])^2d
                  tabdat['lxlbol',j,k] = tabdat['ltotxray',j,k] - $
                     llsun_ergps - tabdat['lbol',j]
                  tabdat['lxlbol_errlo',j,k] = tabdat['ltotxray_errlo',j,k]
                  tabdat['lxlbol_errhi',j,k] = tabdat['ltotxray_errhi',j,k]
               endif
               printf,lun_tmp,gal[i],softstring,hardstring,$
                  format='(A-12,2A27)'
               k++
            endif else begin
               k = 0 ; re-zero x-ray component index
            endelse
         endif
      endif else begin
         k = 0 ; re-zero x-ray component index
      endelse
   endfor

   free_lun,lun_tmp

; Compute physical quantities

;  Median alpha_ox
   igdalphaox = where(tabdat['alphaox'] ne bad)
   print,'Median alpha_ox: ',median(tabdat['alphaox',igdalphaox]),format='(A0,D0.2)'


;  Set AGN fraction to 1 if we don't have a measurement ... but  save  original values first for  table later
   tabdatagnfrac = tabdat['agnfrac']
   tabdatagnfraclb = tabdat['agnfraclb']
   tabdatagnfracub = tabdat['agnfracub']
   ibdagnfrac = where(tabdat['agnfrac'] eq bad,ctbdagnfrac)
   if ctbdagnfrac gt 0 then begin
      tabdat['agnfrac',ibdagnfrac] = 1d
      tabdat['agnfraclb',ibdagnfrac] = 1d
      tabdat['agnfracub',ibdagnfrac] = 1d
   endif

;  Compute L_(F)IR/L_bol
   lirlbol = dblarr(ncos)+bad
   lfirlbol = dblarr(ncos)+bad
   igdlirlbol = where(tabdat['lbol'] ne bad and tabdat['lir'] ne bad)
   igdlfirlbol = where(tabdat['lbol'] ne bad and tabdat['lfir'] ne bad)
   lirlbol[igdlirlbol] = tabdat['lir',igdlirlbol] - tabdat['lbol',igdlirlbol]
   lfirlbol[igdlfirlbol] = tabdat['lfir',igdlfirlbol] - tabdat['lbol',igdlfirlbol]
   
;  Compute L_AGN
   lagn = dblarr(ncos)+bad
   lagnlb = dblarr(ncos)+bad
   lagnub = dblarr(ncos)+bad
   lagn_errlo = dblarr(ncos)+bad
   lagn_errhi = dblarr(ncos)+bad
   igdlbol = where(tabdat['lbol'] ne bad)
   lagn[igdlbol] = tabdat['lbol',igdlbol] + $
                   alog10(tabdat['agnfrac',igdlbol])
   lagnlb[igdlbol] = tabdat['lbol',igdlbol] + $
                     alog10(tabdat['agnfraclb',igdlbol])
   lagnub[igdlbol] = tabdat['lbol',igdlbol] + $
                     alog10(tabdat['agnfracub',igdlbol])
   lagn_errlo[igdlbol] = lagn[igdlbol] - lagnlb[igdlbol]
   lagn_errhi[igdlbol] = lagnub[igdlbol] - lagn[igdlbol]


            
;  Black hole masses
   lmbh = dblarr(ncos)+bad
   lmbh_errlo = dblarr(ncos)+bad
   lmbh_errhi = dblarr(ncos)+bad
   igdrev = where(tabdat['lmbh_rev'] ne bad)
   igdrevrm = where(tabdat['lmbh_rev'] ne bad AND $
                    tabdat['lmbh_rev_type'] eq 'RM')
   igdrevse = where(tabdat['lmbh_rev'] ne bad AND $
                    tabdat['lmbh_rev_type'] eq 'SE')
   igdoth = where(tabdat['lmbh_rev'] ne bad AND $
                  tabdat['lmbh_rev_type'] eq 'GRAVITY')
   ibdrev = where(tabdat['lmbh_rev'] eq bad AND $
                  tabdat['mbh_phot'] ne bad,ctbdrev)
   lmbh[igdrev] = tabdat['lmbh_rev',igdrev]
   lmbh_errlo[igdrevrm] = $
      sqrt(tabdat['lmbh_rev_errlo',igdrevrm]^2d + 0.43d^2d)
   lmbh_errhi[igdrevrm] = $
      sqrt(tabdat['lmbh_rev_errhi',igdrevrm]^2d + 0.43d^2d)
   lmbh_errlo[igdrevse] = $
      sqrt(tabdat['lmbh_rev_errlo',igdrevse]^2d + 2d*0.43d^2d)
   lmbh_errhi[igdrevse] = $
      sqrt(tabdat['lmbh_rev_errhi',igdrevse]^2d + 2d*0.43d^2d)
   lmbh_errlo[igdoth] = -tabdat['lmbh_rev_errlo',igdoth]
   lmbh_errhi[igdoth] = tabdat['lmbh_rev_errhi',igdoth]
   if ctbdrev gt 0 then begin
      lmbh[ibdrev] = alog10(tabdat['mbh_phot',ibdrev])
      lmbh_errlo[ibdrev] = 0.5d
      lmbh_errhi[ibdrev] = 0.5d
   endif
; 
;  Compute Eddington ratio: L_bol/L_Edd = dMdt,acc / dM/dt,Edd
;  L_Edd = eddfac * M_BH
   eddrat = dblarr(ncos)+bad
   leddrat = dblarr(ncos)+bad
   leddrat_errlo = dblarr(ncos)
   leddrat_errhi = dblarr(ncos)
   igd_eddrat = where(lmbh ne bad)
   eddrat[igd_eddrat] = $
      10d^(lagn[igd_eddrat]+llsun_ergps)/$
      (eddfac * 10d^lmbh[igd_eddrat])
   leddrat[igd_eddrat] = alog10(eddrat[igd_eddrat])
   leddrat_errlo[igd_eddrat] = $
      sqrt((lagn[igd_eddrat]-lagnlb[igd_eddrat])^2d + $
           lmbh_errhi[igd_eddrat]^2d)
   leddrat_errhi[igd_eddrat] = $
      sqrt((lagnub[igd_eddrat]-lagn[igd_eddrat])^2d + $
           lmbh_errlo[igd_eddrat]^2d)

;  lamLlam1125
   ldist = lumdist(tabdat['z'],H0=69.3d,Omega_M=0.287d,Lambda0=0.713d,/silent)
   Mpc2cm = IDLUNIT('1d6 pc -> cm')
   fm_unred,1125d*(1d + tabdat['z']),flam1125,tabdat['ebvgal'],flam1125unred
   print,'flam1125 correction for Gal. dust: '
   print,'   median: ',$
      median(flam1125unred/flam1125),format='(A0,D0.2)'
   print,'   max: ',$
      max(flam1125unred/flam1125),format='(A0,D0.2)'
;  factor of 1e-14 is flux unit in table
;  doing computation in observed frame, since flam still in observed frame
   lamLlam1125 = alog10(flam1125unred*4d*!DPi)+2d*alog10(ldist) - 14d + $
      2d*alog10(Mpc2cm.quantity) + alog10(1125d*(tabdat['z']+1d))
   tabdat['luv'] = lamLlam1125
   print,'Median lamLlam1125: ',median(lamLlam1125),format='(A0,D0.2)'

;  Get components
   maxncomp = 15
   nv = orderedhash()
   nv['ncomp'] = intarr(ncos)
   nv['weq'] = dblarr(ncos,3)
   nv['weq_lim'] = dblarr(ncos) + bad ;,3)
   nv['weq_A'] = dblarr(ncos,3)
   nv['weq_lim_A'] = dblarr(ncos) + bad ;,3)
   nv['weq',*,0] = bad
   ;nv['weq_lim',*,0] = bad
   nv['weq_A',*,0] = bad
   ;nv['weq_lim_A',*,0] = bad
   nv['vwtavg'] = dblarr(ncos,3)
   nv['vwtrms'] = dblarr(ncos,3)
   nv['vwtavg',*,0] = bad
   nv['vwtrms',*,0] = bad
   nv['v50'] = dblarr(ncos,maxncomp)+bad
   nv['cf'] = dblarr(ncos,maxncomp)+bad
   nv['tau'] = dblarr(ncos,maxncomp)+bad
   nv['sig'] = dblarr(ncos,maxncomp)+bad
   ovi = orderedhash()
   ovi['ncomp'] = intarr(ncos)
   ovi['weq'] = dblarr(ncos,3)
   ovi['weq_lim'] = dblarr(ncos) + bad ;,3)
   ovi['weq_A'] = dblarr(ncos,3)
   ovi['weq_lim_A'] = dblarr(ncos) + bad ;,3)
   ovi['weq',*,0] = bad
   ;ovi['weq_lim',*,0] = bad
   ovi['weq_A',*,0] = bad
   ;ovi['weq_lim_A',*,0] = bad
   ovi['vwtavg'] = dblarr(ncos,3)
   ovi['vwtrms'] = dblarr(ncos,3)
   ovi['vwtavg',*,0] = bad
   ovi['vwtrms',*,0] = bad
   ovi['v50'] = dblarr(ncos,maxncomp)+bad
   ovi['cf'] = dblarr(ncos,maxncomp)+bad
   ovi['tau'] = dblarr(ncos,maxncomp)+bad
   ovi['sig'] = dblarr(ncos,maxncomp)+bad
   pv = orderedhash()
   pv['ncomp'] = intarr(ncos)
   pv['weq'] = dblarr(ncos,3)
   pv['weq_A'] = dblarr(ncos,3)
   pv['weq',*,0] = bad
   pv['weq_A',*,0] = bad
   pv['vwtavg'] = dblarr(ncos,3)
   pv['vwtrms'] = dblarr(ncos,3)
   pv['vwtavg',*,0] = bad
   pv['vwtrms',*,0] = bad
   pv['v50'] = dblarr(ncos,maxncomp)+bad
   pv['cf'] = dblarr(ncos,maxncomp)+bad
   pv['tau'] = dblarr(ncos,maxncomp)+bad
   pv['sig'] = dblarr(ncos,maxncomp)+bad
   for i=0,ncos-1 do begin
      fitdir_gal=fitdir+tabdat['sgal',i]+'/'
      if file_test(fitdir_gal,/dir) then begin
         file_tmp = fitdir_gal+tabdat['sgal',i]+'NVpar_best.txt'
         xdr_tmp = fitdir_gal+tabdat['sgal',i]+'NV_fit.xdr'
         if file_test(file_tmp) then begin
            readcol,file_tmp,ncomp_tmp,/silent,numline=1,skip=3,format='(I0,X)'
;            readcol,file_tmp,ncompem_tmp,/silent,numline=1,skip=1,format='(I0,X,X,X,X)'
            readcol,file_tmp,totweq_tmp,totweq_errlo_tmp,totweq_errhi_tmp,$
               /silent,numline=1,skip=4,format='(D0,D0,D0,X)'
            readcol,file_tmp,vwtavg_tmp,vwtavg_errlo_tmp,vwtavg_errhi_tmp,$
               /silent,numline=1,skip=5,format='(D0,D0,D0,X)'
            readcol,file_tmp,vwtrms_tmp,vwtrms_errlo_tmp,vwtrms_errhi_tmp,$
               /silent,numline=1,skip=6,format='(D0,D0,D0,X)'
            nv['ncomp',i]=ncomp_tmp[0]
            if totweq_tmp[0] gt 0 then begin
               nv['weq',i,0]=alog10(totweq_tmp[0])
               nv['weq',i,1]=alog10(totweq_tmp[0])-$
                  alog10(totweq_tmp[0]-totweq_errlo_tmp[0])
               nv['weq',i,2]=alog10(totweq_tmp[0]+totweq_errhi_tmp[0])-$
                  alog10(totweq_tmp[0])
               nv['weq_A',i,0]=totweq_tmp[0]
               nv['weq_A',i,1]=totweq_errlo_tmp[0]
               nv['weq_A',i,2]=totweq_errhi_tmp[0]
            endif
            nv['vwtavg',i,0] = vwtavg_tmp[0]
            nv['vwtavg',i,1] = vwtavg_errlo_tmp[0]
            nv['vwtavg',i,2] = vwtavg_errhi_tmp[0]
            nv['vwtrms',i,0] = vwtrms_tmp[0]
            nv['vwtrms',i,1] = vwtrms_errlo_tmp[0]
            nv['vwtrms',i,2] = vwtrms_errhi_tmp[0]
            readcol,file_tmp,cf,tau1243,lambda1243,sig,vel,/silent,$
                    numline=ncomp_tmp[0],skip=10,format='(D,D,D,D,D)'
            nv['v50',i,0:ncomp_tmp[0]-1] = vel
            nv['cf',i,0:ncomp_tmp[0]-1] = cf
            nv['tau',i,0:ncomp_tmp[0]-1] = tau1243
            nv['sig',i,0:ncomp_tmp[0]-1] = sig
         endif
         file_tmp = fitdir_gal+tabdat['sgal',i]+'OVIpar_best.txt'
         if file_test(file_tmp) then begin
            readcol,file_tmp,ncomp_tmp,/silent,numline=1,skip=3,format='(I0,X)'
;            readcol,file_tmp,ncompem_tmp,/silent,numline=1,skip=1,format='(I0,X,X,X,X)'
            readcol,file_tmp,totweq_tmp,totweq_errlo_tmp,totweq_errhi_tmp,$
               /silent,numline=1,skip=4,format='(D0,D0,D0,X)'
            readcol,file_tmp,vwtavg_tmp,vwtavg_errlo_tmp,vwtavg_errhi_tmp,$
               /silent,numline=1,skip=5,format='(D0,D0,D0,X)'
            readcol,file_tmp,vwtrms_tmp,vwtrms_errlo_tmp,vwtrms_errhi_tmp,$
               /silent,numline=1,skip=6,format='(D0,D0,D0,X)'
            ovi['ncomp',i]=ncomp_tmp[0]
            if totweq_tmp[0] gt 0 then begin
               ovi['weq',i,0]=alog10(totweq_tmp[0])
               ovi['weq',i,1]=alog10(totweq_tmp[0])-$
                  alog10(totweq_tmp[0]-totweq_errlo_tmp[0])
               ovi['weq',i,2]=alog10(totweq_tmp[0]+totweq_errhi_tmp[0])-$
                  alog10(totweq_tmp[0])
               ovi['weq_A',i,0]=totweq_tmp[0]
               ovi['weq_A',i,1]=totweq_errlo_tmp[0]
               ovi['weq_A',i,2]=totweq_errhi_tmp[0]
            endif
            ovi['vwtavg',i,0] = vwtavg_tmp[0]
            ovi['vwtavg',i,1] = vwtavg_errlo_tmp[0]
            ovi['vwtavg',i,2] = vwtavg_errhi_tmp[0]
            ovi['vwtrms',i,0] = vwtrms_tmp[0]
            ovi['vwtrms',i,1] = vwtrms_errlo_tmp[0]
            ovi['vwtrms',i,2] = vwtrms_errhi_tmp[0]
            readcol,file_tmp,cf,tau1038,lambda1038,sig,vel,/silent,$
                    numline=ncomp_tmp[0],skip=10,format='(D,D,D,D,D)'
            ovi['v50',i,0:ncomp_tmp[0]-1] = vel
            ovi['cf',i,0:ncomp_tmp[0]-1] = cf
            ovi['tau',i,0:ncomp_tmp[0]-1] = tau1038
            ovi['sig',i,0:ncomp_tmp[0]-1] = sig
         endif
         file_tmp = fitdir_gal+tabdat['sgal',i]+'PVpar_best.txt'
         if file_test(file_tmp) then begin
            readcol,file_tmp,ncomp_tmp,/silent,numline=1,skip=3,format='(I0,X)'
            ;            readcol,file_tmp,ncompem_tmp,/silent,numline=1,skip=1,format='(I0,X,X,X,X)'
            readcol,file_tmp,totweq_tmp,totweq_errlo_tmp,totweq_errhi_tmp,$
               /silent,numline=1,skip=4,format='(D0,D0,D0,X)'
            readcol,file_tmp,vwtavg_tmp,vwtavg_errlo_tmp,vwtavg_errhi_tmp,$
               /silent,numline=1,skip=5,format='(D0,D0,D0,X)'
            readcol,file_tmp,vwtrms_tmp,vwtrms_errlo_tmp,vwtrms_errhi_tmp,$
               /silent,numline=1,skip=6,format='(D0,D0,D0,X)'
            pv['ncomp',i]=ncomp_tmp[0]
            if totweq_tmp[0] gt 0 then begin
               pv['weq',i,0]=alog10(totweq_tmp[0])
               pv['weq',i,1]=alog10(totweq_tmp[0])-$
                  alog10(totweq_tmp[0]-totweq_errlo_tmp[0])
               pv['weq',i,2]=alog10(totweq_tmp[0]+totweq_errhi_tmp[0])-$
                  alog10(totweq_tmp[0])
               pv['weq_A',i,0]=totweq_tmp[0]
               pv['weq_A',i,1]=totweq_errlo_tmp[0]
               pv['weq_A',i,2]=totweq_errhi_tmp[0]
            endif
            pv['vwtavg',i,0] = vwtavg_tmp[0]
            pv['vwtavg',i,1] = vwtavg_errlo_tmp[0]
            pv['vwtavg',i,2] = vwtavg_errhi_tmp[0]
            pv['vwtrms',i,0] = vwtrms_tmp[0]
            pv['vwtrms',i,1] = vwtrms_errlo_tmp[0]
            pv['vwtrms',i,2] = vwtrms_errhi_tmp[0]
            readcol,file_tmp,cf,tau1128,lambda1128,sig,vel,/silent,$
               numline=ncomp_tmp[0],skip=10,format='(D,D,D,D,D)'
            pv['v50',i,0:ncomp_tmp[0]-1] = vel
            pv['cf',i,0:ncomp_tmp[0]-1] = cf
            pv['tau',i,0:ncomp_tmp[0]-1] = tau1128
            pv['sig',i,0:ncomp_tmp[0]-1] = sig
         endif
      endif
      ; compute upper limits
      iwin_A = 0.5d ; half-width of window in which to compute stats
      ispec_A = 2.5d ; half-width of window in which to compute spectrum
      modsig = 50d ; sigma for upper limit model
      if nv['ncomp',i] eq 0 OR ovi['ncomp',i] eq 0 then begin
         readcol,specdir+tabdat['sgal',i]+'.txt',lam,fx,err,format='(D,D,D)',$
            /silent
         npts = n_elements(lam)
         disp = median(lam[1:npts-1]-lam[0:npts-2])
         iwin_p = round(iwin_A / disp)
         ispec_p = round(ispec_A / disp)
         lovi1032 = (1d + tabdat['z',i])*linelist['OVI1031']
         lovi1038 = (1d + tabdat['z',i])*linelist['OVI1037']
         lnv1238 = (1d + tabdat['z',i])*linelist['NV1238']
         lnv1242 = (1d + tabdat['z',i])*linelist['NV1242']
         ; to avoid strong Galactic absorption for this case
         if tabdat['sgal',i] eq 'pg0157' then begin
            lovi1032 -= 2d
            lovi1038 -= 2d
         endif
         iovi1032 = value_locate(lam,lovi1032)
         iovi1038 = value_locate(lam,lovi1038)
         inv1238 = value_locate(lam,lnv1238)
         inv1242 = value_locate(lam,lnv1242)
         ; note we remove PG1501 b/c no data in NV region
         if nv['ncomp',i] eq 0 AND inv1238 ne -1 AND inv1238 ne npts-1 AND $
            inv1242 ne -1 AND inv1242 ne npts-1 AND $
            tabdat['sgal',i] ne 'pg1501'then begin
            medfx = median([fx[inv1238-iwin_p:inv1238+iwin_p],$
               fx[inv1242-iwin_p:inv1242+iwin_p]])
            mederr = median([err[inv1238-iwin_p:inv1238+iwin_p],$
               err[inv1242-iwin_p:inv1242+iwin_p]])
            rmsfx = stddev([fx[inv1238-iwin_p:inv1238+iwin_p],$
               fx[inv1242-iwin_p:inv1242+iwin_p]])
            mederrnorm = mederr / medfx
            rmsfxnorm = rmsfx / medfx
            nvmodweq = 1b
            ;nvmodweqlo = 1b
            ;nvmodweqhi = 1b
            ; model is: cf = rmsfxnorm/2 (divide by 2 b/c two lines!)
            ;    tau = 5
            ;    rest wave
            ;    median measured sig
            nvmod = ifsf_doubletfcn(lam[inv1238-ispec_p:inv1242+ispec_p],$
               [1,0,rmsfxnorm/2d,5d,lnv1242,modsig],doubletname='NV',$
               weq=nvmodweq)
            ;nvmodlo = ifsf_doubletfcn(lam[inv1238-ispec_p:inv1242+ispec_p],$
            ;   [1,0,rmsfxnorm/4d,5d,lnv1242,modsig],doubletname='NV',$
            ;   weq=nvmodweqlo)
            ;nvmodhi = ifsf_doubletfcn(lam[inv1238-ispec_p:inv1242+ispec_p],$
            ;   [1,0,rmsfxnorm,5d,lnv1242,modsig],doubletname='NV',$
            ;   weq=nvmodweqhi)
            nv['weq_lim',i] = alog10(nvmodweq.abs[0])
            ;nv['weq_lim',i,0] = alog10(nvmodweq.abs[0])
            ;nv['weq_lim',i,1] = alog10(nvmodweq.abs[0]) - $
            ;   alog10(nvmodweqlo.abs[0])
            ;nv['weq_lim',i,2] = alog10(nvmodweqhi.abs[0]) - $
            ;   alog10(nvmodweq.abs[0])
            nv['weq_lim_A',i] = nvmodweq.abs[0]
            ;nv['weq_lim_A',i,0] = nvmodweq.abs[0]
            ;nv['weq_lim_A',i,1] = nvmodweq.abs[0] - nvmodweqlo.abs[0]
            ;nv['weq_lim_A',i,2] = nvmodweqhi.abs[0] - nvmodweq.abs[0]
            ; Examine model spectra cf actual spectra
            ;set_plot,'x'
            ;cgplot,lam[inv1238-ispec_p:inv1242+ispec_p],nvmod,yran=[0,1.5]
            ;cgoplot,lam[inv1238-ispec_p:inv1242+ispec_p],fx[inv1238-ispec_p:inv1242+ispec_p]/medfx
            ;print,tabdat['sgal',i],mederrnorm,nvmodweq.abs[0]
            
         endif
         ; note we remove PG2349 b/c geocoronal Lyalpha
         if ovi['ncomp',i] eq 0 AND iovi1032 ne -1 AND iovi1032 ne npts-1 AND $
            iovi1038 ne -1 AND iovi1038 ne npts-1 AND $
            tabdat['sgal',i] ne 'pg2349' then begin
            medfx = median([fx[iovi1032-iwin_p:iovi1032+iwin_p],$
               fx[iovi1038-iwin_p:iovi1038+iwin_p]])
            mederr = median([err[iovi1032-iwin_p:iovi1032+iwin_p],$
               err[iovi1038-iwin_p:iovi1038+iwin_p]])
            rmsfx = stddev([fx[iovi1032-iwin_p:iovi1032+iwin_p],$
               fx[iovi1038-iwin_p:iovi1038+iwin_p]])
            mederrnorm = mederr / medfx
            rmsfxnorm = rmsfx / medfx
            ovimodweq = 1b
            ;ovimodweqlo = 1b
            ;ovimodweqhi = 1b
            ; model is: cf = rmsfxnorm/2 (divide by 2 b/c two lines!)
            ;    tau = 5
            ;    rest wave
            ;    median measured sig
            ovimod = ifsf_doubletfcn(lam[iovi1032-ispec_p:iovi1038+ispec_p],$
               [1,0,rmsfxnorm/2d,5d,lovi1038,modsig],doubletname='OVI',$
               weq=ovimodweq)
            ;ovimodlo = ifsf_doubletfcn(lam[iovi1032-ispec_p:iovi1038+ispec_p],$
            ;   [1,0,rmsfxnorm/4d,5d,lovi1038,modsig],doubletname='OVI',$
            ;   weq=ovimodweqlo)
            ;ovimodhi = ifsf_doubletfcn(lam[iovi1032-ispec_p:iovi1038+ispec_p],$
            ;   [1,0,rmsfxnorm,5d,lovi1038,modsig],doubletname='OVI',$
            ;   weq=ovimodweqhi)
            ovi['weq_lim',i] = alog10(ovimodweq.abs[0])
            ;ovi['weq_lim',i,0] = alog10(ovimodweq.abs[0])
            ;ovi['weq_lim',i,1] = alog10(ovimodweq.abs[0]) - $
            ;   alog10(ovimodweqlo.abs[0])
            ;ovi['weq_lim',i,2] = alog10(ovimodweqhi.abs[0]) - $
            ;   alog10(ovimodweq.abs[0])
            ovi['weq_lim_A',i] = ovimodweq.abs[0]
            ;ovi['weq_lim_A',i,0] = ovimodweq.abs[0]
            ;ovi['weq_lim_A',i,1] = ovimodweq.abs[0] - ovimodweqlo.abs[0]
            ;ovi['weq_lim_A',i,2] = ovimodweqhi.abs[0] - ovimodweq.abs[0]
            ; Examine model spectra cf actual spectra
            ;set_plot,'x'
            ;cgplot,lam[iovi1032-ispec_p:iovi1038+ispec_p],ovimod,yran=[0,1.5]
            ;cgoplot,lam[iovi1032-ispec_p:iovi1038+ispec_p],fx[iovi1032-ispec_p:iovi1038+ispec_p]/medfx
            ;print,tabdat['sgal',i],mederrnorm,ovimodweq.abs[0]
         endif
     endif
   endfor
   igd_nv_comp = where(nv['v50'] ne bad,ctnvcomp)
   igd_ovi_comp = where(ovi['v50'] ne bad,ctovicomp)

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d


;  Plots

; Colors
tvlct,[[27],[158],[119]],100
tvlct,[[217],[95],[2]],101
tvlct,[[117],[112],[179]],102
tvlct,[[231],[41],[138]],103
tvlct,[[102],[166],[30]],104
tvlct,[[230],[171],[2]],105
tvlct,[[166],[118],[29]],106
tvlct,[[102],[102],[102]],107
tvlct,[[0],[0],[0]],108

; Arrays for non-detection plots

inv_nd = where(tabdat['nvstatus'] eq 'U',ctnv_nd)
iovi_nd = where(tabdat['ovistatus'] eq 'U',ctovi_nd)
tmpy_nv = dblarr(ctnv_nd)+1d99
tmpy_ovi = dblarr(ctovi_nd)+1d99
for i=0,ctnv_nd-1 do tmpy_nv[i] *= randomu(seed)+0.5d
for i=0,ctovi_nd-1 do tmpy_ovi[i] *= randomu(seed)+0.5d

;  output table for regressions
; Threshold pval for marking as significant in output regression table
tpval = 0.05
openw,lun_stat,tabdir+'tab_regressions.txt',/get_lun
printf,lun_stat,'#Col 1: x-axis quantity'
printf,lun_stat,'#Col 2: y-axis quantity'
printf,lun_stat,'#Col 3-4: p, lower limit flag'
printf,lun_stat,'#Col 6-8: r, -/+ error'
printf,lun_stat,'#Col 9: no. of points'
statform = '(A10,A10,D7.3,I2,3D6.2,I3)'

openw,lun_reg_tex,tabdir+'tab_regressions.tex',/get_lun
tweq = '$W_{\rm eq}$'
tvel = ['$v_{\rm wtavg}$','$\sigma_{\rm rms}$']
regtexform = '(A20,A3,A70,A3,I3,A3,D7.3,A3,'+$
   'D6.2,A0,D-4.2,A0,D-4.2,A0,A3)'
bregtexform = '(A20,A3,A11,A58,A1,A3,I3,A3,D7.3,A3,'+$
      'D6.2,A0,D-4.2,A0,D-4.2,A0,A3)'
regtexform_pvallim = '(A20,A3,A70,A3,I3,A3,A3,D-5.3,A2,'+$
   'D6.2,A0,D-4.2,A0,D-4.2,A0,A3)'
bregtexform_pvallim = '(A20,A3,A11,A58,A1,A3,I3,A3,A3,D-5.3,A2,'+$
      'D6.2,A0,D-4.2,A0,D-4.2,A0,A3)'
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'z'
   xtit = 'Redshift'
   xran = [0d,0.399d]

   ylab = 'lbol'
   ytit = 'log(L$\downbol$/L$\sun$)'
   yran = [11.2,13.2]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=2,default=4,/quiet
   cgplot,tabdat[xlab],tabdat[ylab],xran=xran,yran=yran,$
          psym=9,symsize=1.5,color='Black',xtit=xtit,ytit=ytit,position=[0.18,0.18,0.99,0.99]

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; v50 vs.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.2,13.2]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-8000,1000]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=9,symsize=1.5,$
           color='Blue'
;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'luv'
   xtit = 'log[ $\lambda$L$\down\\lambda$ (1125 $\Angstrom$)/erg s$\up-1$]'
   xran = [42.5d,47d]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-8000,1000]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.02]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,position=[0.18,0.18,0.99,0.99]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.2,13.2]

   xrebin = rebin(lagn,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,position=[0.18,0.18,0.99,0.99]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.499]

   xrebin = rebin(lmbh,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,position=[0.18,0.18,0.99,0.99]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.499]

   xrebin = rebin(leddrat,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,position=[0.18,0.18,0.99,0.99]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'alphaox'
   xtit = '$\alpha$$\downox$'
   xran = [-2.199,-1.1]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,position=[0.18,0.18,0.99,0.99]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=9,symsize=1.5,$
           color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close
   
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'cf'
   xtit = 'Covering Factor'
   xran = [-0.05,1.05]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,nv[xlab,igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          position=[0.18,0.18,0.99,0.99]
;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,ovi[xlab,igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=9,symsize=1.5,$
           color='Blue'
;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; sig vs
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.2,13.2]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
   yran = [0.5,3.5]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)
   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          position=[0.14,0.14,0.95,0.95]
;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=9,symsize=1.5,$
           color='Blue'
;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'luv'
   xtit = 'log[ $\lambda$L$\down\\lambda$ (1125 $\Angstrom$)/erg s$\up-1$]'
   xran = [42.5d,47d]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
      yran = [0.5,3.5]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)
   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
      position=[0.14,0.14,0.95,0.95]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.05]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)
   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          position=[0.14,0.14,0.95,0.95]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.3,12.8]

   xrebin = rebin(lagn,ncos,maxncomp)
   xrebin = rebin(lagn,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          position=[0.14,0.14,0.95,0.95]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.5]

   xrebin = rebin(lmbh,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          position=[0.14,0.14,0.95,0.95]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.5]

   xrebin = rebin(leddrat,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          position=[0.14,0.14,0.95,0.95]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=9,symsize=1.5,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close
   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'alphaox'
   xtit = '$\alpha$$\downox$'
   xran = [-2.199,-1.1]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          position=[0.14,0.14,0.95,0.95]
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=9,symsize=1.5,$
           color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close
   
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'cf'
   xtit = 'Covering Factor'
   xran = [-0.1,1.1]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,nv[xlab,igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          position=[0.14,0.14,0.95,0.95]
;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,ovi[xlab,igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=9,symsize=1.5,$
           color='Blue'
;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; weq vs.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   chars = 2.5

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xtex = 'log($L_{\rm BOL}/L_\sun$)'
   xran = [11.2,13.2]

   ylab = 'weq'
   ytit = 'log(W$\downeq$ / $\Angstrom$)'
   yran = [-2,2]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,tabdat[xlab],nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d
   cgoplot,tabdat[xlab],ovi['weq_lim'],symsize=1,color='Blue',$
           psym='Open Down Triangle'


   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      tabdat[xlab] ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      tabdat[xlab] ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      tabdat[xlab] ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      tabdat[xlab] ne bad,ctovilim)

   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = tabdat[xlab,iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'luv'
   xtit = 'log[ $\lambda$L$\down\\lambda$ (1125 $\Angstrom$)/erg s$\up-1$]'
   xtex = 'log[$\lambda L_{1125}/{\rm erg~s}^{-1}$]'
   xran = [42.5d,47d]


   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,tabdat[xlab],nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d
   cgoplot,tabdat[xlab],ovi['weq_lim'],symsize=1,color='Blue',$
           psym='Open Down Triangle'


   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      tabdat[xlab] ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      tabdat[xlab] ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      tabdat[xlab] ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      tabdat[xlab] ne bad,ctovilim)

   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = tabdat[xlab,iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   ;xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]


   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,yerr=yerr,/metro,$ ;xerr=xerr,
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         ;xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xtex = xtit
   xran = [0.65,1.02]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab,*,0],xran=xran,yran=yran,$
      psym=1,symsize=0.01,color='Grey',xtit=xtit,$
      position=[0.15,0.15,0.95,0.95],$
      err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d,$
      err_xlo=tabdat[xlab]-tabdat[xlab+'lb'],$
      err_xhi=tabdat[xlab+'ub'] - tabdat[xlab],/err_clip,$
      xtickint=0.1
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,tabdat[xlab],nv['weq_lim'],symsize=0.01,color='Grey',psym=1,$
      err_xlo=tabdat[xlab]-tabdat[xlab+'lb'],$
      err_xhi=tabdat[xlab+'ub'] - tabdat[xlab],/err_clip,$
      err_width=0d
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=1,symsize=0.01,color='Grey',$
      err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d,$
      err_xlo=tabdat[xlab]-tabdat[xlab+'lb'],$
      err_xhi=tabdat[xlab+'ub'] - tabdat[xlab],/err_clip
   cgoplot,tabdat[xlab],ovi['weq_lim'],symsize=0.01,color='Grey',psym=1,$
      err_xlo=tabdat[xlab]-tabdat[xlab+'lb'],$
      err_xhi=tabdat[xlab+'ub'] - tabdat[xlab],/err_clip,$
      err_width=0d
   cgoplot,tabdat[xlab],nv[ylab,*,0],psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab],nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=9,symsize=1.5,color='Blue'
   cgoplot,tabdat[xlab],ovi['weq_lim'],symsize=1,color='Blue',$
      psym='Open Down Triangle'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      tabdat[xlab] ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      tabdat[xlab] ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      tabdat[xlab] ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      tabdat[xlab] ne bad,ctovilim)

   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = tabdat[xlab,iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   xerravg = (-tabdat[xlab+'lb']+tabdat[xlab+'ub'])/2d
   xerr = xerravg[iall]
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d] ;mean(yerravg[invdet])]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=3,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.35*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.15,0.40,0.40,0.15,0.15],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.15,0.7,0.4,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xtex = 'log($L_{\rm AGN}/L_\sun$)'
   xran = [11.301,12.799]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,lagn,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=1,symsize=0.01,color='Grey',xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d,$
          err_xlo=lagn_errlo,$
          err_xhi=lagn_errhi,/err_clip,err_color='Grey',$
          xtickint=0.4
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,lagn,nv['weq_lim'],symsize=0.01,color='Grey',psym=1,$
      err_xlo=lagn_errlo,$
      err_xhi=lagn_errhi,/err_clip,$
      err_width=0d,err_color='Grey'
   cgoplot,lagn,ovi[ylab,*,0],psym=1,symsize=0.01,color='Blue',$
      err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d,$
      err_xlo=lagn_errlo,$
      err_xhi=lagn_errhi,/err_clip,err_color='Grey'
   cgoplot,lagn,ovi['weq_lim'],symsize=0.01,color='Blue',psym=1,$
      err_xlo=lagn_errlo,$
      err_xhi=lagn_errhi,/err_clip,$
      err_width=0d,err_color='Grey'
   cgoplot,lagn,nv[ylab,*,0],psym=15,symsize=1,color='Red'
   cgoplot,lagn,nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,lagn,ovi[ylab,*,0],psym=9,symsize=1.5,color='Blue'
   cgoplot,lagn,ovi['weq_lim'],symsize=1,color='Blue',$
      psym='Open Down Triangle'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      lagn ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      lagn ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      lagn ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      lagn ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      lagn ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      lagn ne bad,ctovilim)
   
   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = lagn[iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   xerravg = (lagn_errlo+lagn_errhi)/2d
   xerr = xerravg[iall]
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d] ;mean(yerravg[invdet])]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xtex=  'log($M_{\rm BH}/M_\sun$)'
   xran = [6.5,9.5]


   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,lmbh,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=1,symsize=0.01,color='Grey',xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d,$
          err_xlo=lmbh_errlo,$
          err_xhi=lmbh_errhi,/err_clip,err_color='Grey',$
          xtickint=0.4
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,lmbh,nv['weq_lim'],symsize=0.01,color='Grey',psym=1,$
      err_xlo=lmbh_errlo,$
      err_xhi=lmbh_errhi,/err_clip,$
      err_width=0d,err_color='Grey'
   cgoplot,lmbh,ovi[ylab,*,0],psym=1,symsize=0.01,color='Blue',$
      err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d,$
      err_xlo=lmbh_errlo,$
      err_xhi=lmbh_errhi,/err_clip,err_color='Grey'
   cgoplot,lmbh,ovi['weq_lim'],symsize=0.01,color='Blue',psym=1,$
      err_xlo=lmbh_errlo,$
      err_xhi=lmbh_errhi,/err_clip,$
      err_width=0d,err_color='Grey'
   cgoplot,lmbh,nv[ylab,*,0],psym=15,symsize=1,color='Red'
   cgoplot,lmbh,nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,lmbh,ovi[ylab,*,0],psym=9,symsize=1.5,color='Blue'
   cgoplot,lmbh,ovi['weq_lim'],symsize=1,color='Blue',$
      psym='Open Down Triangle'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      lmbh ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      lmbh ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      lmbh ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      lmbh ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      lmbh ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      lmbh ne bad,ctovilim)
   
   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = lmbh[iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   xerravg = (lmbh_errlo+lmbh_errhi)/2d
   xerr = xerravg[iall]
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d] ;mean(yerravg[invdet])]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.35*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.15,0.40,0.40,0.15,0.15],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.15,0.7,0.40,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xtex = xtit
   xran = [-2,0.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,leddrat,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=1,symsize=0.01,color='Grey',xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d,$
          err_xlo=leddrat_errlo,$
          err_xhi=leddrat_errhi,/err_clip,err_color='Grey',$
          xtickint=0.4
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,leddrat,nv['weq_lim'],symsize=0.01,color='Grey',psym=1,$
      err_xlo=leddrat_errlo,$
      err_xhi=leddrat_errhi,/err_clip,$
      err_width=0d,err_color='Grey'
   cgoplot,leddrat,ovi[ylab,*,0],psym=1,symsize=0.01,color='Blue',$
      err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d,$
      err_xlo=leddrat_errlo,$
      err_xhi=leddrat_errhi,/err_clip,err_color='Grey'
   cgoplot,leddrat,ovi['weq_lim'],symsize=0.01,color='Blue',psym=1,$
      err_xlo=leddrat_errlo,$
      err_xhi=leddrat_errhi,/err_clip,$
      err_width=0d,err_color='Grey'
   cgoplot,leddrat,nv[ylab,*,0],psym=15,symsize=1,color='Red'
   cgoplot,leddrat,nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,leddrat,ovi[ylab,*,0],psym=9,symsize=1.5,color='Blue'
   cgoplot,leddrat,ovi['weq_lim'],symsize=1,color='Blue',$
      psym='Open Down Triangle'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      leddrat ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      leddrat ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      leddrat ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      leddrat ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      leddrat ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      leddrat ne bad,ctovilim)
   
   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = leddrat[iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   xerravg = (leddrat_errlo+leddrat_errhi)/2d
   xerr = xerravg[iall]
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d] ;mean(yerravg[invdet])]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'alphaox'
   xtit = '$\alpha$$\downox$'
   xtex = '$\alpha_{\rm OX}$
   xran = [-2.2,-1.1]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,tabdat[xlab],nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d
   cgoplot,tabdat[xlab],ovi['weq_lim'],symsize=1,color='Blue',$
           psym='Open Down Triangle'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      tabdat[xlab] ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      tabdat[xlab] ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      tabdat[xlab] ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      tabdat[xlab] ne bad,ctovilim)

   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = tabdat[xlab,iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lirlbol'
   xtit = 'log(L$\downIR$/L$\downbol$)'
   xtex = 'log($L_{\rm IR}/L_{\rm BOL}$)'
   xran = [-0.799,0.199]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,lirlbol,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,lirlbol,nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,lirlbol,ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d
   cgoplot,lirlbol,ovi['weq_lim'],symsize=1,color='Blue',$
           psym='Open Down Triangle'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      lirlbol ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      lirlbol ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      lirlbol ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      lirlbol ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      lirlbol ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      lirlbol ne bad,ctovilim)

   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = lirlbol[iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.35*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.15,0.4,0.4,0.15,0.15],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.15,0.7,0.4,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lfirlbol'
   xtit = 'log(L$\downFIR$/L$\downbol$)'
   xtex = 'log($L_{\rm FIR}/L_{\rm BOL}$)'
   xran = [-2,0]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,lfirlbol,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,lfirlbol,nv['weq_lim'],symsize=1,color='Red',$
      psym='Filled Down Triangle'
   cgoplot,lfirlbol,ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d
   cgoplot,lfirlbol,ovi['weq_lim'],symsize=1,color='Blue',$
           psym='Open Down Triangle'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      lfirlbol ne bad,ctbothdet)
   ibothlim = where(nv['weq_lim'] ne bad AND $
      ovi['weq_lim'] ne bad AND $
      lfirlbol ne bad,ctbothlim)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      lfirlbol ne bad,ctnvdet)
   invlim = where(nv['weq_lim'] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      ovi['weq_lim'] eq bad AND $
      lfirlbol ne bad,ctnvlim)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      lfirlbol ne bad,ctovidet)
   iovilim = where(ovi['weq_lim'] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      nv['weq_lim'] eq bad AND $
      lfirlbol ne bad,ctovilim)

   npts = ctbothdet + ctbothlim + ctnvdet + ctnvlim + ctovidet + ctovilim
   iall = [ibothdet,ibothlim,invdet,invlim,iovidet,iovilim]
   xdat = lfirlbol[iall]
   ydat = [alog10((10d^nv[ylab,ibothdet,0]+10d^ovi[ylab,ibothdet,0])/2d),$
      alog10((10d^nv['weq_lim',ibothlim]+10d^ovi['weq_lim',ibothlim])/2d),$
      nv[ylab,invdet,0],nv['weq_lim',invlim],$
      ovi[ylab,iovidet,0],ovi['weq_lim',iovilim]]
   xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],dblarr(ctbothlim)+0.1d,$
      ynverravg[invdet],dblarr(ctnvlim)+0.1d,$
      yovierravg[iovidet],dblarr(ctovilim)+0.1d]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctbothlim),$
      bytarr(ctnvdet)+1b,bytarr(ctnvlim),$
      bytarr(ctovidet)+1b,bytarr(ctovilim)]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'nhxray'
   xtit = 'log[ N(H, X-ray) / cm$\up-2$]'
   xtex = 'log[$N({\rm H})/{\rm cm}^{-2}$]'
   xran = [19.5d,24.5d]

   openw,tmplun,plotdir+ylab+'_vs_'+xlab+'_dat.txt',/get_lun
   printf,tmplun,'xdat','ydat','xerr','yerr','xl','yl',$
      format='(4A8,2A3)'

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         xerrlo = tabdat[xlab+'_errlo',i,0:tabdat['n_xray',i]-1]
         xerrhi = tabdat[xlab+'_errhi',i,0:tabdat['n_xray',i]-1]
         x_lim = tabdat[xlab+'_lim',i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i,0]],tabdat['n_xray',i])
         ynv_errlo = rebin([nv[ylab,i,1]],tabdat['n_xray',i])
         ynv_errhi = rebin([nv[ylab,i,2]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i,0]],tabdat['n_xray',i])
         yovi_errlo = rebin([ovi[ylab,i,1]],tabdat['n_xray',i])
         yovi_errhi = rebin([ovi[ylab,i,2]],tabdat['n_xray',i])
         ynv_lim = rebin([nv[ylab+'_lim',i]],tabdat['n_xray',i])
         yovi_lim = rebin([ovi[ylab+'_lim',i]],tabdat['n_xray',i])

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x_lim,ynv,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Left Triangle',/err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         usersym,[-0.366d,1d,0d,-0.366d]*1.5d,[-0.366d,0d,1d,-0.366d]*1.5d,$
            thick=!P.thick,/fill
         cgoplot,x_lim,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym=8,/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x_lim,yovi,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Left Triangle',/err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         usersym,[-0.366d,1d,0d,-0.366d]*1.5d,[-0.366d,0d,1d,-0.366d]*1.5d,$
            thick=!P.thick
         cgoplot,x_lim,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym=8,/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'
                  
         ; this assumes all X-ray obs. for an object are either detections or limits
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
            xuselim = 0b
         endif else if x_lim[0] ne bad then begin
            xuse = mean(x_lim[0])
            xuseerr = 0.3d ; assume factor-of-2 error on limit
            xuselim = 1b
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 0b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 0b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 1b
            endif else if yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 1b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_ovi OR not nodat_nv then begin
            printf,tmplun,xuse,yuse,xuseerr,yuseerr,xuselim,yuselim,$
            format='(D8.4,D8.4,D8.4,D8.4,I3,I3)'
            xdat = [xdat, xuse]
            ydat = [ydat, yuse]
            cens = [cens, xuselim OR yuselim]
         endif
      endif
   endfor

   readcol,plotdir+ylab+'_vs_'+xlab+'_stat.txt',cc,pval,pvaldist,skip=1,/silent,$
      format='(D,D,A)'
   ; translate pymcc results to LINMIX_ERR results
   pval[0] = double(pvaldist[0])
   if pvaldist[1] eq 'True' then begin
      pval[1] = 1b
      pvalstr = 'p<'
      pvallimstr = '$<$'
   endif else begin
      pval[1] = 0b
      pvalstr = 'p='
      pvallimstr = ''
   endelse
   corr = dblarr(3)
   corr[0] = cc[1]
   corr[1] = cc[1]-cc[0]
   corr[2] = cc[2]-cc[1]
   xloc = xran[0]+0.35*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string(pvalstr,pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 0b)
   indet = where(cens eq 1b)
   cgpolygon,[0.15,0.4,0.4,0.15,0.15],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.15,0.7,0.4,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),$
         amp,pvallimstr,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform_pvallim
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pvallimstr,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform_pvallim
   endelse

   cgps_close

   free_lun,tmplun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'gamxray'
   xtit = '$\Gamma$ (X-ray)'
   xtex = '$\Gamma$'
   xran = [1,3.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         xerrlo = tabdat[xlab+'_errlo',i,0:tabdat['n_xray',i]-1]
         xerrhi = tabdat[xlab+'_errhi',i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i,0]],tabdat['n_xray',i])
         ynv_errlo = rebin([nv[ylab,i,1]],tabdat['n_xray',i])
         ynv_errhi = rebin([nv[ylab,i,2]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i,0]],tabdat['n_xray',i])
         yovi_errlo = rebin([ovi[ylab,i,1]],tabdat['n_xray',i])
         yovi_errhi = rebin([ovi[ylab,i,2]],tabdat['n_xray',i])
         ynv_lim = rebin([nv[ylab+'_lim',i]],tabdat['n_xray',i])
         yovi_lim = rebin([ovi[ylab+'_lim',i]],tabdat['n_xray',i])

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 0b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR not yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 0b
            endif else if not yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 0b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif
   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fsoftxray'
   xtit = 'log[F(0.5-2 keV)/erg s$\up-1$ cm$\up-2$]'
   xtex = 'log[$F(0.5-2~{\rm keV})/{\rm erg~s}^{-1}~{\rm cm}^{-2}$]'
   xran = [-14,-10]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      ; special treatment for F(soft); PG0050 has one obs. with no flux listed
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)
         ynv_lim = rebin([nv[ylab+'_lim',i]],ngdx)
         yovi_lim = rebin([ovi[ylab+'_lim',i]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 0b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR not yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 0b
            endif else if not yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 0b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif
   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fhardxray'
   xtit = 'log[F(2-10 keV)/erg s$\up-1$ cm$\up-2$]'
   xtex = 'log[$F(2-10~{\rm keV})/{\rm erg~s}^{-1}~{\rm cm}^{-2}$]'
   xran = [-13.3,-9.8]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],xtickint=1d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      ; special treatment for F(soft); PG0050 has one obs. with no flux listed
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)
         ynv_lim = rebin([nv[ylab+'_lim',i]],ngdx)
         yovi_lim = rebin([ovi[ylab+'_lim',i]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 0b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR not yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 0b
            endif else if not yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 0b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif
   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftxray'
   xtit = 'log[L(0.5-2 keV)/erg s$\up-1$]'
   xtex = 'log[$L(0.5-2~{\rm keV})/{\rm erg~s}^{-1}$]'
   xran = [42,46]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],xtickint=1d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      ; special treatment for F(soft); PG0050 has one obs. with no flux listed
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)
         ynv_lim = rebin([nv[ylab+'_lim',i]],ngdx)
         yovi_lim = rebin([ovi[ylab+'_lim',i]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 0b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR not yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 0b
            endif else if not yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 0b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif
   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lhardxray'
   xtit = 'log[L(2-10 keV)/erg s$\up-1$]'
   xtex = 'log[$L(2-10~{\rm keV})/{\rm erg~s}^{-1}$]'
   xran = [42.25,46.25]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],xtickint=1d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      ; special treatment for F(soft); PG0050 has one obs. with no flux listed
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)
         ynv_lim = rebin([nv[ylab+'_lim',i]],ngdx)
         yovi_lim = rebin([ovi[ylab+'_lim',i]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 0b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR not yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 0b
            endif else if not yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 0b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif
   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'ltotxray'
   xtit = 'log[L(0.5-10 keV)/erg s$\up-1$]'
   xtex = 'log[$L(0.5-10~{\rm keV})/{\rm erg~s}^{-1}$]'
   xran = [42.5,46.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95],xtickint=1d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      ; special treatment for F(soft); PG0050 has one obs. with no flux listed
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)
         ynv_lim = rebin([nv[ylab+'_lim',i]],ngdx)
         yovi_lim = rebin([ovi[ylab+'_lim',i]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 0b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR not yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 0b
            endif else if not yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 0b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif
   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftratxray'
   xtit = 'L(0.5-2 keV) / L(0.5-10 keV)'
   xtex = 'log[$L(0.5-2~{\rm keV})/L(0.5-10~{\rm keV})$]'
   xran = [-0.8d,0d]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      ; special treatment for F(soft); PG0050 has one obs. with no flux listed
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_err',i,igdx]
         xerrhi = tabdat[xlab+'_err',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)
         ynv_lim = rebin([nv[ylab+'_lim',i]],ngdx)
         yovi_lim = rebin([ovi[ylab+'_lim',i]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 0b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR not yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 0b
            endif else if not yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 0b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif
   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lxlbol'
   xtit = 'log[L(0.5-10 keV) / L$\downbol$]'
   xtex = 'log[$L(0.5-10~{\rm keV})/L_{\rm BOL}$]'
   xran = [-3.3,-0.3]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.15,0.15,0.95,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      ; special treatment for F(soft); PG0050 has one obs. with no flux listed
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)
         ynv_lim = rebin([nv[ylab+'_lim',i]],ngdx)
         yovi_lim = rebin([ovi[ylab+'_lim',i]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else if ynv_lim[0] ne bad then begin
            yuse = ynv_lim[0]
            yuseerr = 0.1d
            yuselim = 0b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv OR not yuselim then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (yovi[0] + ynv[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else if yovi_lim[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi_lim[0]
               yuseerr = 0.1d
               yuselim = 0b
            endif else if not yuselim then begin
               yuse = (yovi_lim[0] + ynv_lim[0])/2d
               yuseerr = 0.1d
               yuselim = 0b
            endif
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv_lim,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Down Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi_lim,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Down Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif
   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.7,0.95,0.95,0.7,0.7],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.7,0.7,0.95,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; vwtavg OR vwtrms vs.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   chars=2.5

   vlabs = ['vwtavg','vwtrms']
   vtits = ['v$\down wtavg$ (km/s)','$\sigma$$\down rms$ (km/s)']
   vrans = [[1000,-7000],[-100,3500]]

   for vind=0,1 do begin

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xtex = 'log($L_{\rm BOL}/L_\sun$)'
   xran = [11.2,13.2]

   ylab = vlabs[vind]
   ytit = vtits[vind]
   yran = vrans[*,vind]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.25,0.15,0.99,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      tabdat[xlab] ne bad,ctbothdet)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctnvdet)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctovidet)

   npts = ctbothdet + ctnvdet + ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = tabdat[xlab,iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   xlab = 'luv'
   xtit = 'log[ $\lambda$L$\down\\lambda$ (1125 $\Angstrom$)/erg s$\up-1$]'
   xtex = 'log[$\lambda L_{1125}/{\rm erg~s}^{-1}$]'
   xran = [42.5d,46.99d]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.25,0.15,0.99,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      tabdat[xlab] ne bad,ctbothdet)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctnvdet)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctovidet)

   npts = ctbothdet + ctnvdet + ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = tabdat[xlab,iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   ;xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,yerr=yerr,/metro,$ ;xerr=xerr,
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         ;xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'alphaox'
   xtit = '$\alpha$$\downox$'
   xtex = '$\alpha_{\rm OX}$
   xran = [-2.199,-1.101]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.25,0.15,0.99,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      tabdat[xlab] ne bad,ctbothdet)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctnvdet)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      tabdat[xlab] ne bad,ctovidet)

   npts = ctbothdet + ctnvdet + ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = tabdat[xlab,iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   ;xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,yerr=yerr,/metro,$ ;xerr=xerr,
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         ;xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lirlbol'
   xtit = 'log(L$\downIR$/L$\downbol$)'
   xtex = 'log($L_{\rm IR}/L_{\rm BOL}$)'
   xran = [-0.799,0.199]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,lirlbol,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.25,0.15,0.99,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,lirlbol,ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      lirlbol ne bad,ctbothdet)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      lirlbol ne bad,ctnvdet)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      lirlbol ne bad,ctovidet)

   npts = ctbothdet + ctnvdet + ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = lirlbol[iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lfirlbol'
   xtit = 'log(L$\downFIR$/L$\downbol$)'
   xtex = 'log($L_{\rm FIR}/L_{\rm BOL}$)'
   xran = [-2,0]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,lfirlbol,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,$
          position=[0.25,0.15,0.99,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,lfirlbol,ovi[ylab,*,0],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      lfirlbol ne bad,ctbothdet)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      lfirlbol ne bad,ctnvdet)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      lfirlbol ne bad,ctovidet)

   npts = ctbothdet + ctnvdet + ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = lfirlbol[iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   xerr = dblarr(npts)+0.1d
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close



   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xtex = xtit
   xran = [0.65,1.02]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab,*,0],xran=xran,yran=yran,$
      psym=1,symsize=0.01,color='Grey',xtit=xtit,$
      position=[0.25,0.15,0.99,0.95],$
      err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d,$
      err_xlo=tabdat[xlab]-tabdat[xlab+'lb'],$
      err_xhi=tabdat[xlab+'ub'] - tabdat[xlab],/err_clip,$
      xtickint=0.1d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=1,symsize=0.01,color='Grey',$
      err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d,$
      err_xlo=tabdat[xlab]-tabdat[xlab+'lb'],$
      err_xhi=tabdat[xlab+'ub'] - tabdat[xlab],/err_clip
   cgoplot,tabdat[xlab],nv[ylab,*,0],psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab],ovi[ylab,*,0],psym=9,symsize=1.5,color='Blue'

   ; NV + OVI stats
   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND tabdat[xlab] ne bad,ctbothdet)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND tabdat[xlab] ne bad,ctnvdet)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND tabdat[xlab] ne bad,ctovidet)
   npts = ctbothdet+ctnvdet +ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = tabdat[xlab,iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   xerravg = (-tabdat[xlab+'lb']+tabdat[xlab+'ub'])/2d
   xerr = xerravg[iall]
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   redolinmix = 0b
   ; special case; ngauss = 3 won't converge quickly for vwtrms
   if vind eq 0 then ngauss = 3 else ngauss = 1
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=ngauss,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xtex = 'log($L_{\rm AGN}/L_\sun$)'
   xran = [11.301,12.799]
   
   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,lagn,nv[ylab,*,0],xran=xran,yran=yran,$
      psym=1,symsize=0.01,color='Grey',xtit=xtit,$
      position=[0.25,0.15,0.99,0.95],$ ;,xtickform='(A1)',$
      err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d,$
      err_xlo=lagn_errlo,$
      err_xhi=lagn_errhi,/err_clip,err_color='Grey',xtickint=0.4d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,lagn,ovi[ylab,*,0],psym=1,symsize=0.01,color='Blue',$
      err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d,$
      err_xlo=lagn_errlo,$
      err_xhi=lagn_errhi,/err_clip,err_color='Grey'
   cgoplot,lagn,nv[ylab,*,0],psym=15,symsize=1,color='Red'
   cgoplot,lagn,ovi[ylab,*,0],psym=9,symsize=1.5,color='Blue'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      lagn ne bad,ctbothdet)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      lagn ne bad,ctnvdet)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      lagn ne bad,ctovidet)

   npts = ctbothdet + ctnvdet + ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = lagn[iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   xerravg = (lagn_errlo+lagn_errhi)/2d
   xerr = xerravg[iall]
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xtex=  'log($M_{\rm BH}/M_\sun$)'
   xran = [6.501,9.499]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,lmbh,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=1,symsize=0.01,color='Grey',xtit=xtit,$
          position=[0.25,0.15,0.99,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d,$
          err_xlo=lmbh_errlo,$
          err_xhi=lmbh_errhi,/err_clip,err_color='Grey'
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,lmbh,ovi[ylab,*,0],psym=1,symsize=0.01,color='Blue',$
      err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d,$
      err_xlo=lmbh_errlo,$
      err_xhi=lmbh_errhi,/err_clip,err_color='Grey'
   cgoplot,lmbh,nv[ylab,*,0],psym=15,symsize=1,color='Red'
   cgoplot,lmbh,ovi[ylab,*,0],psym=9,symsize=1.5,color='Blue'

   ibothdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] ne bad AND $
      lmbh ne bad,ctbothdet)
   invdet = where(nv[ylab,*,0] ne bad AND $
      ovi[ylab,*,0] eq bad AND $
      lmbh ne bad,ctnvdet)
   iovidet = where(ovi[ylab,*,0] ne bad AND $
      nv[ylab,*,0] eq bad AND $
      lmbh ne bad,ctovidet)

   npts = ctbothdet + ctnvdet + ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = lmbh[iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   xerravg = (lmbh_errlo+lmbh_errhi)/2d
   xerr = xerravg[iall]
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xtex = xtit
   xran = [-1.999,0.499]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,leddrat,nv[ylab,*,0],xran=xran,yran=yran,$
          psym=1,symsize=0.01,color='Grey',xtit=xtit,$
          position=[0.25,0.15,0.99,0.95],$ ;,xtickform='(A1)',$
          err_ylo=nv[ylab,*,1],err_yhi=nv[ylab,*,2],err_width=0d,$
          err_xlo=leddrat_errlo,$
          err_xhi=leddrat_errhi,/err_clip,err_color='Grey'
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   cgoplot,leddrat,ovi[ylab,*,0],psym=1,symsize=0.01,color='Blue',$
      err_ylo=ovi[ylab,*,1],err_yhi=ovi[ylab,*,2],err_width=0d,$
      err_xlo=leddrat_errlo,$
      err_xhi=leddrat_errhi,/err_clip,err_color='Grey'
   cgoplot,leddrat,nv[ylab,*,0],psym=15,symsize=1,color='Red'
   cgoplot,leddrat,ovi[ylab,*,0],psym=9,symsize=1.5,color='Blue'

   npts = ctbothdet + ctnvdet + ctovidet
   iall = [ibothdet,invdet,iovidet]
   xdat = leddrat[iall]
   ydat = [(nv[ylab,ibothdet,0]+ovi[ylab,ibothdet,0])/2d,$
      nv[ylab,invdet,0],ovi[ylab,iovidet,0]]
   xerravg = (leddrat_errlo+leddrat_errhi)/2d
   xerr = xerravg[iall]
   ynverravg = (nv[ylab,*,1] + nv[ylab,*,2])/2d
   yovierravg = (ovi[ylab,*,1] + ovi[ylab,*,2])/2d
   ybotherravg = (ynverravg + yovierravg)/2d
   yerr = [ybotherravg[ibothdet],ynverravg[invdet],yovierravg[iovidet]]
   cens = [bytarr(ctbothdet)+1b,bytarr(ctnvdet)+1b,bytarr(ctovidet)+1b]

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 1b)
   indet = where(cens eq 0b)
   cgpolygon,[0.74,0.99,0.99,0.74,0.74],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.74,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   ;cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'nhxray'
   xtit = 'log[ N(H, X-ray) / cm$\up-2$]'
   xtex = 'log[$N({\rm H})/{\rm cm}^{-2}$]'
   xran = [19.5d,24.5d]

   openw,tmplun,plotdir+ylab+'_vs_'+xlab+'_dat.txt',/get_lun
   printf,tmplun,'xdat','ydat','xerr','yerr','xl','yl',$
      format='(4A8,2A3)'

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
      position=[0.25,0.15,0.99,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         xerrlo = tabdat[xlab+'_errlo',i,0:tabdat['n_xray',i]-1]
         xerrhi = tabdat[xlab+'_errhi',i,0:tabdat['n_xray',i]-1]
         x_lim = tabdat[xlab+'_lim',i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i,0]],tabdat['n_xray',i])
         ynv_errlo = rebin([nv[ylab,i,1]],tabdat['n_xray',i])
         ynv_errhi = rebin([nv[ylab,i,2]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i,0]],tabdat['n_xray',i])
         yovi_errlo = rebin([ovi[ylab,i,1]],tabdat['n_xray',i])
         yovi_errhi = rebin([ovi[ylab,i,2]],tabdat['n_xray',i])

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x_lim,ynv,symsize=1.5,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            psym='Filled Left Triangle',/err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x_lim,yovi,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            psym='Open Left Triangle',/err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

         ; this assumes all X-ray obs. for an object are either detections or limits
         nodat_nv = 0b
         nodat_ovi = 0b
         yuselim = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
            xuselim = 0b
         endif else if x_lim[0] ne bad then begin
            xuse = mean(x_lim[0])
            xuseerr = 0.3d ; assume factor-of-2 error on limit
            xuselim = 1b
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad AND nodat_nv then begin
            yuse = yovi[0]
            yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_ovi OR not nodat_nv then begin
            printf,tmplun,xuse,yuse,xuseerr,yuseerr,xuselim,yuselim,$
               format='(D8.4,D12.2,D8.4,D7.3,I3,I3)'
            xdat = [xdat, xuse]
            ydat = [ydat, yuse]
            cens = [cens, xuselim]
         endif
      endif
   endfor

   readcol,plotdir+ylab+'_vs_'+xlab+'_stat.txt',cc,pval,pvaldist,skip=1,/silent,$
      format='(D,D,A)'
   ; translate pymcc results to LINMIX_ERR results
   pval[0] = double(pvaldist[0])
   if pvaldist[1] eq 'True' then begin
      pval[1] = 1b
      pvalstr = 'p<'
      pvallimstr = '$<$'
   endif else begin
      pval[1] = 0b
      pvalstr = 'p='
      pvallimstr = ''
   endelse
   corr = dblarr(3)
   corr[0] = cc[1]
   corr[1] = cc[1]-cc[0]
   corr[2] = cc[2]-cc[1]
   xloc = xran[0]+0.35*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string(pvalstr,pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   idet = where(cens eq 0b)
   indet = where(cens eq 1b)
   cgpolygon,[0.15,0.4,0.4,0.15,0.15],[0.7,0.7,0.95,0.95,0.7],$
      /fill,fcol='white',/norm
   cgplot,xdat[idet],ydat[idet],psym=16,symsize=0.75,color='Black',/noerase,$
      pos=[0.15,0.7,0.4,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      xticks=1,yticks=1,xminor=1,yminor=1
   cgoplot,xdat[indet],ydat[indet],psym=9,symsize=0.75,color='Black'


   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),$
         amp,pvallimstr,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform_pvallim
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tweq,amp,xtex,amp,n_elements(xdat),amp,pvallimstr,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform_pvallim
   endelse

   cgps_close

   free_lun,tmplun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   xlab = 'gamxray'
   xtit = '$\Gamma$ (X-ray)'
   xtex = '$\Gamma$'
   xran = [1,3.499]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.25,0.15,0.99,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (ynv[0] + yovi[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif

   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   cgplot,xdat,ydat,psym=16,symsize=1,color='Black',/noerase,$
      pos=[0.75,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      background='White'

   cgps_close

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fsoftxray'
   xtit = 'log[F(0.5-2 keV)/erg s$\up-1$ cm$\up-2$]'
   xtex = 'log[$F(0.5-2~{\rm keV})/{\rm erg~s}^{-1}~{\rm cm}^{-2}$]'
   xran = [-13.99,-10.001]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.25,0.15,0.99,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (ynv[0] + yovi[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif

   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   cgplot,xdat,ydat,psym=16,symsize=1,color='Black',/noerase,$
      pos=[0.75,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      background='White'

   cgps_close

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fhardxray'
   xtit = 'log[F(2-10 keV)/erg s$\up-1$ cm$\up-2$]'
   xtex = 'log[$F(2-10~{\rm keV})/{\rm erg~s}^{-1}~{\rm cm}^{-2}$]'
   xran = [-13.3,-9.8]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.25,0.15,0.99,0.95],xtickint=1d
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (ynv[0] + yovi[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif

   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   cgplot,xdat,ydat,psym=16,symsize=1,color='Black',/noerase,$
      pos=[0.75,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      background='White'

   cgps_close

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftxray'
   xtit = 'log[L(0.5-2 keV)/erg s$\up-1$]'
   xtex = 'log[$L(0.5-2~{\rm keV})/{\rm erg~s}^{-1}$]'
   xran = [42,45.999]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.25,0.15,0.99,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (ynv[0] + yovi[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif

   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   cgplot,xdat,ydat,psym=16,symsize=1,color='Black',/noerase,$
      pos=[0.75,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      background='White'

   cgps_close

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lhardxray'
   xtit = 'log[L(2-10 keV)/erg s$\up-1$]'
   xtex = 'log[$L(2-10~{\rm keV})/{\rm erg~s}^{-1}$]'
   xran = [42.25,46.25]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.25,0.15,0.99,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (ynv[0] + yovi[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif

   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1] - 0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   cgplot,xdat,ydat,psym=16,symsize=1,color='Black',/noerase,$
      pos=[0.75,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      background='White'

   cgps_close

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'ltotxray'
   xtit = 'log[L(0.5-10 keV)/erg s$\up-1$]'
   xtex = 'log[$L(0.5-10~{\rm keV})/{\rm erg~s}^{-1}$]'
   xran = [42.5,46.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.25,0.15,0.99,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (ynv[0] + yovi[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif

   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   cgplot,xdat,ydat,psym=16,symsize=1,color='Black',/noerase,$
      pos=[0.75,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      background='White'

   cgps_close

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftratxray'
   xtit = 'L(0.5-2 keV) / L(0.5-10 keV)'
   xtex = 'log[$L(0.5-2~{\rm keV})/L(0.5-10~{\rm keV})$]'
   xran = [-0.8d,0d]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.25,0.15,0.99,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_err',i,igdx]
         xerrhi = tabdat[xlab+'_err',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (ynv[0] + yovi[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif

   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   cgplot,xdat,ydat,psym=16,symsize=1,color='Black',/noerase,$
      pos=[0.75,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      background='White'

   cgps_close

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lxlbol'
   xtit = 'log[L(0.5-10 keV) / L$\downbol$]'
   xtex = 'log[$L(0.5-10~{\rm keV})/L_{\rm BOL}$]'
   xran = [-3.3,-0.3]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,$
          position=[0.25,0.15,0.99,0.95]
   cgtext,ytit,0.07,0.55,/norm,align=0.5,orient=90d
   xdat = !NULL
   ydat = !NULL
   xerr = !NULL
   yerr = !NULL
   cens = !NULL
   for i=0,ncos-1 do begin
      igdx = where(tabdat[xlab,i,0:tabdat['n_xray',i]-1] ne bad,ngdx)
      if ngdx gt 0 then begin
         x = tabdat[xlab,i,igdx]
         xerrlo = tabdat[xlab+'_errlo',i,igdx]
         xerrhi = tabdat[xlab+'_errhi',i,igdx]
         ynv = rebin([nv[ylab,i,0]],ngdx)
         ynv_errlo = rebin([nv[ylab,i,1]],ngdx)
         ynv_errhi = rebin([nv[ylab,i,2]],ngdx)
         yovi = rebin([ovi[ylab,i,0]],ngdx)
         yovi_errlo = rebin([ovi[ylab,i,1]],ngdx)
         yovi_errhi = rebin([ovi[ylab,i,2]],ngdx)

         ; this assumes all X-ray obs. for an object are detections
         nodat_nv = 0b
         nodat_ovi = 0b
         if x[0] ne bad then begin
            xuse = mean(x)
            xuseerr = (mean(xerrlo) + mean(xerrhi)) / 2d
         endif else begin
            nodat_nv=1b
            nodat_ovi=1b
         endelse
         if ynv[0] ne bad then begin
            yuse = ynv[0]
            yuseerr = (ynv_errlo[0]+ynv_errhi[0])/2d
            yuselim = 1b
         endif else begin
            nodat_nv = 1b
         endelse
         if yovi[0] ne bad then begin
            if nodat_nv then begin
               yuse = yovi[0]
               yuseerr = (yovi_errlo[0]+yovi_errhi[0])/2d
               yuselim = 1b
            endif else begin
               yuse = (ynv[0] + yovi[0])/2d
               yuseerr = (yuseerr + (yovi_errlo[0]+yovi_errhi[0])/2d)/2d
            endelse
         endif else begin
            nodat_ovi = 1b
         endelse
         if not nodat_nv OR not nodat_ovi then begin
            xdat = [xdat,xuse]
            xerr = [xerr,xuseerr]
            ydat = [ydat,yuse]
            yerr = [yerr,yuseerr]
            cens = [cens,yuselim]
         endif

         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=ynv_errlo,err_yhi=ynv_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_xlo=xerrlo,$
            err_xhi=xerrhi,err_ylo=yovi_errlo,err_yhi=yovi_errhi,err_width=0d,$
            /err_clip
         cgoplot,x,yovi,/linesty,color='Blue'

      endif

   endfor

   ; NV + OVI stats
   redolinmix = 0b
   if redolinmix then begin
      fitout = drt_runlinmix(xdat,ydat,xerr=xerr,yerr=yerr,/metro,$
         alpha=alpha,beta=beta,corr=corr,sigsqr=sigsqr,pval=pval,$
         miniter=lm_miniter,maxiter=lm_maxiter,detected=cens,ngauss=1,/verbose)
      linmixpar = {alpha:alpha,$
         beta:beta,$
         corr:corr,$
         sigsqr:sigsqr,$
         nfit:npts,$
         pval:pval,$
         xdat:xdat,$
         ydat:ydat,$
         xerr:xerr,$
         yerr:yerr,$
         yfit:fitout[*,0],$
         yfitlo:fitout[*,1],$
         yfithi:fitout[*,2],$
         yfit2lo:fitout[*,3],$
         yfit2hi:fitout[*,4]}
      save,linmixpar,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
   endif else begin
      restore,file=plotdir+ylab+'_vs_'+xlab+'_NVOVI.xdr'
      pval = linmixpar.pval
      corr = linmixpar.corr
   endelse

   xloc = xran[0]+0.05*(xran[1]-xran[0])
   yloc = yran[1]-0.05*(yran[1]-yran[0])
   lab = string('p=',pval[0],'  r=',corr[0],$
      format='(A0,D0.3,A0,D0.2)')+$
      textoidl('_{-'+string(corr[1],format='(D0.2)')+'}^{+'+$
      string(corr[2],format='(D0.2)')+'}')+$
      string('  N=',n_elements(xdat),format='(A0,I0)')
   cgtext,xloc,yloc,lab,/dat,chars=1.25

   cgplot,xdat,ydat,psym=16,symsize=1,color='Black',/noerase,$
      pos=[0.75,0.7,0.99,0.95],xran=xran,yran=yran,xtickf='(A1)',ytickf='(A1)',$
      background='White'

   cgps_close

   if pval[0] lt tpval then begin
      printf,lun_stat,xlab.ToUpper(),ylab.ToUpper(),pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform
      printf,lun_reg_tex,tvel[vind],amp,'\underline{ ',xtex,'}',amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=bregtexform
   endif else begin
      printf,lun_stat,xlab,ylab,pval[0],pval[1],corr[0],corr[1],corr[2],$
         n_elements(xdat),format=statform  
      printf,lun_reg_tex,tvel[vind],amp,xtex,amp,n_elements(xdat),amp,pval[0],amp,$
         corr[0],lerr0,corr[1],lerr1,corr[2],lerr2,dslash,format=regtexform
   endelse

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; vwtrms vs. vwtavg
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'vwtavg'
   xtit = 'Average depth-weighted velocity (km/s)'
   xran = [1000,-9000]

   ylab = 'vwtrms'
   ytit = 'Average depth-weighted velocity RMS (km/s)'
   yran = [0,3500]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=chars,default=4,/quiet
   cgplot,nv[xlab,0],nv[ylab,0],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,xtit=xtit,$
          position=[0.18,0.18,0.98,0.98]
   cgoplot,ovi[xlab,0],ovi[ylab,0],psym=9,symsize=1.5,$
           color='Blue'
   ;cgoplot,xran,yran,/linesty
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'v50'
   xtit = 'v$\down50$ (km/s)'
   xran = [1000,-8000]

   ylab = 'quasar_zsort'
   ytit = 'Quasars, sorted by z'
   yran = [0,ncos+1]

;  sort quasars by z
   isort = sort(tabdat['z'])
   gal_sort = tabdat['gal',isort]
   nvstatus_sort = tabdat['nvstatus',isort]
   ovistatus_sort = tabdat['ovistatus',isort]
   nv_v50_sort = nv['v50']
   ovi_v50_sort = ovi['v50']
   nv_sig_sort = nv['sig']
   ovi_sig_sort = ovi['sig']
   for i=0,ncos-1 do begin
      nv_v50_sort[i,*] = nv['v50',isort[i],*]
      ovi_v50_sort[i,*] = ovi['v50',isort[i],*]
      nv_sig_sort[i,*] = nv['sig',isort[i],*]
      ovi_sig_sort[i,*] = ovi['sig',isort[i],*]
   endfor
   yarr = indgen(ncos)+1
   yrebin = rebin(indgen(ncos)+1,ncos,maxncomp)

   ix_nv = where(nvstatus_sort eq 'X',ct_x_nv)
   ix_ovi = where(ovistatus_sort eq 'X',ct_x_ovi)
   iu_nv = where(nvstatus_sort eq 'U',ct_u_nv)
   iu_ovi = where(ovistatus_sort eq 'U',ct_u_ovi)
   ig_nv = where(nvstatus_sort eq 'G' OR nvstatus_sort eq 'GF',ct_g_nv)
   ig_ovi = where(ovistatus_sort eq 'G' OR ovistatus_sort eq 'GF',ct_g_ovi)
   
   igd_nv_comp_sort = where(nv_v50_sort ne bad)
   igd_ovi_comp_sort = where(ovi_v50_sort ne bad)
   sigthresh = 25d
   igd_nv_comp_sort_nar = where(nv_sig_sort ne bad AND nv_sig_sort le sigthresh)
   igd_ovi_comp_sort_nar = where(ovi_sig_sort ne bad AND ovi_sig_sort le sigthresh)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=10,/nomatch,$
             charsize=1,default=4,/quiet

   cgplot,nv_v50_sort[igd_nv_comp_sort],yrebin[igd_nv_comp_sort],$
          xran=xran,yran=yran,xtit=xtit,ytit=ytit,$
          psym=15,symsize=1,color='Red',$
          ytickname=gal_sort,ytickv=indgen(ncos)+1,yticks=ncos-1,$
          position=[0.2,0.15,0.8,0.99],/noerase
   cgoplot,nv_v50_sort[igd_nv_comp_sort_nar],yrebin[igd_nv_comp_sort_nar],$
          psym=6,symsize=2,color='Red'
   cgoplot,ovi_v50_sort[igd_ovi_comp_sort],yrebin[igd_ovi_comp_sort],$
           psym=9,symsize=1.5,color='Blue'
   cgoplot,ovi_v50_sort[igd_ovi_comp_sort_nar],yrebin[igd_ovi_comp_sort_nar],$
           psym=9,symsize=2,color='Blue'
   j=1
   for i=1,ncos do begin
      if j mod 5 eq 0 then cgoplot,xran,[i,i],linesty=1
      j++
   endfor
   
   xran = [0,3]
   cgplot,[0],/nodata,xran=xran,yran=yran,/xsty,/ysty,$
          position=[0.85,0.15,0.99,0.99],/noerase,$
          ytickf='(A1)',xticks=3,xtickv=[1,2],xtickname=['NV','OVI']
   cgoplot,intarr(ct_x_nv)+1,yarr[ix_nv],psym=7,color='Red'
   cgoplot,intarr(ct_x_ovi)+2,yarr[ix_ovi],psym=7,color='Blue'
   cgoplot,intarr(ct_u_nv)+1,yarr[iu_nv],psym=11,color='Red'
   cgoplot,intarr(ct_u_ovi)+2,yarr[iu_ovi],psym=11,color='Blue'
   cgoplot,intarr(ct_g_nv)+1,yarr[ig_nv],psym=36,color='Red',symsize=1.5
   cgoplot,intarr(ct_g_ovi)+2,yarr[ig_ovi],psym=36,color='Blue',symsize=1.5
 
   j=1
   for i=1,ncos do begin
      if j mod 5 eq 0 then cgoplot,xran,[i,i],linesty=1
      j++
   endfor

   !P.position=[0,1,0,1]

;   al_legend,['NV 1238,1243','O VI 1032,1038','not in spectral range',$
;      'affected by geocoronal line or chip gap','undetected'],$
;      /norm,spacing=1.5,$
;      color=['Red','Blue','Black','Black','Black'],$
;      psym=[6,16,7,36,11],symsize=[1.5,1,1,1.5,1],position=[0.555,0.12]
   al_legend,['NV 1238,1243','O VI 1032,1038',$
              'NV, $\sigma$ < 25 km/s',$
              'O VI, $\sigma$ < 25 km/s'],$
              /norm,spacing=1.5,$
              color=['Red','Blue','Red','Blue'],$
              psym=[15,9,6,9],symsize=[1,1.5,2,2],position=[0.2,0.1],margin=.6
   al_legend,['not in spectral range',$
              'affected by geocoronal line or chip gap','undetected'],$
              /norm,spacing=1.5,$
              color=['Black','Black','Black'],$
              psym=[7,36,11],symsize=[1,1.5,1],position=[0.5,0.1],margin=.6

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'v50'
   xtit = 'v$\down50$ (km/s)'
   xran = [1000,-8000]

   ylab = 'quasar_lbolsort'
   ytit = 'Quasars, sorted by L$\downbol$'
   yran = [0,ncos+1]

;  sort quasars by z
   isort = sort(tabdat['lbol'])
   gal_sort = tabdat['gal',isort]
   nvstatus_sort = tabdat['nvstatus',isort]
   ovistatus_sort = tabdat['ovistatus',isort]
   nv_v50_sort = nv['v50']
   ovi_v50_sort = ovi['v50']
   nv_sig_sort = nv['sig']
   ovi_sig_sort = ovi['sig']
   for i=0,ncos-1 do begin
      nv_v50_sort[i,*] = nv['v50',isort[i],*]
      ovi_v50_sort[i,*] = ovi['v50',isort[i],*]
      nv_sig_sort[i,*] = nv['sig',isort[i],*]
      ovi_sig_sort[i,*] = ovi['sig',isort[i],*]
   endfor
   yarr = indgen(ncos)+1
   yrebin = rebin(indgen(ncos)+1,ncos,maxncomp)

   ix_nv = where(nvstatus_sort eq 'X',ct_x_nv)
   ix_ovi = where(ovistatus_sort eq 'X',ct_x_ovi)
   iu_nv = where(nvstatus_sort eq 'U',ct_u_nv)
   iu_ovi = where(ovistatus_sort eq 'U',ct_u_ovi)
   ig_nv = where(nvstatus_sort eq 'G' OR nvstatus_sort eq 'GF',ct_g_nv)
   ig_ovi = where(ovistatus_sort eq 'G' OR ovistatus_sort eq 'GF',ct_g_ovi)
   
   igd_nv_comp_sort = where(nv_v50_sort ne bad)
   igd_ovi_comp_sort = where(ovi_v50_sort ne bad)
   sigthresh = 25d
   igd_nv_comp_sort_nar = where(nv_sig_sort ne bad AND nv_sig_sort le sigthresh)
   igd_ovi_comp_sort_nar = where(ovi_sig_sort ne bad AND ovi_sig_sort le sigthresh)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=10,/nomatch,$
             charsize=1,default=4,/quiet

   cgplot,nv_v50_sort[igd_nv_comp_sort],yrebin[igd_nv_comp_sort],$
          xran=xran,yran=yran,xtit=xtit,ytit=ytit,$
          psym=15,symsize=1,color='Red',$
          ytickname=gal_sort,ytickv=indgen(ncos)+1,yticks=ncos-1,$
          position=[0.2,0.15,0.8,0.99],/noerase
   cgoplot,nv_v50_sort[igd_nv_comp_sort_nar],yrebin[igd_nv_comp_sort_nar],$
          psym=6,symsize=2,color='Red'
   cgoplot,ovi_v50_sort[igd_ovi_comp_sort],yrebin[igd_ovi_comp_sort],$
           psym=9,symsize=1.5,color='Blue'
   cgoplot,ovi_v50_sort[igd_ovi_comp_sort_nar],yrebin[igd_ovi_comp_sort_nar],$
           psym=9,symsize=2,color='Blue'
   j=1
   for i=1,ncos do begin
      if j mod 5 eq 0 then cgoplot,xran,[i,i],linesty=1
      j++
   endfor
   
   xran = [0,3]
   cgplot,[0],/nodata,xran=xran,yran=yran,/xsty,/ysty,$
          position=[0.85,0.15,0.99,0.99],/noerase,$
          ytickf='(A1)',xticks=3,xtickv=[1,2],xtickname=['NV','OVI']
   cgoplot,intarr(ct_x_nv)+1,yarr[ix_nv],psym=7,color='Red'
   cgoplot,intarr(ct_x_ovi)+2,yarr[ix_ovi],psym=7,color='Blue'
   cgoplot,intarr(ct_u_nv)+1,yarr[iu_nv],psym=11,color='Red'
   cgoplot,intarr(ct_u_ovi)+2,yarr[iu_ovi],psym=11,color='Blue'
   cgoplot,intarr(ct_g_nv)+1,yarr[ig_nv],psym=36,color='Red',symsize=1.5
   cgoplot,intarr(ct_g_ovi)+2,yarr[ig_ovi],psym=36,color='Blue',symsize=1.5
 
   j=1
   for i=1,ncos do begin
      if j mod 5 eq 0 then cgoplot,xran,[i,i],linesty=1
      j++
   endfor

   !P.position=[0,1,0,1]

;   al_legend,['NV 1238,1243','O VI 1032,1038','not in spectral range',$
;      'affected by geocoronal line or chip gap','undetected'],$
;      /norm,spacing=1.5,$
;      color=['Red','Blue','Black','Black','Black'],$
;      psym=[6,16,7,36,11],symsize=[1.5,1,1,1.5,1],position=[0.555,0.12]
   al_legend,['NV 1238,1243','O VI 1032,1038',$
              'NV, $\sigma$ < 25 km/s',$
              'O VI, $\sigma$ < 25 km/s'],$
              /norm,spacing=1.5,$
              color=['Red','Blue','Red','Blue'],$
              psym=[15,9,6,9],symsize=[1,1.5,2,2],position=[0.2,0.1],margin=.6
   al_legend,['not in spectral range',$
              'affected by geocoronal line or chip gap','undetected'],$
              /norm,spacing=1.5,$
              color=['Black','Black','Black'],$
              psym=[7,36,11],symsize=[1,1.5,1],position=[0.5,0.1],margin=.6

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'hist_cf.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.5,default=4,/quiet

   bin1 = cghistogram(nv['cf',igd_nv_comp],binsize=0.2d,min=0d,missing=bad,/freq)
   bin2 = cghistogram(ovi['cf',igd_ovi_comp],binsize=0.2d,min=0d,missing=bad,/freq)
   bin3 = (bin1+bin2)/2d
;   move 1s into top bin ...
   bin1[4]+=bin1[5]
   bin2[4]+=bin2[5]
   bin3[4]+=bin3[5]
   bin1n = cghistogram(nv['cf',igd_nv_comp],binsize=0.2d,min=0d,max=1d,missing=bad)
   bin2n = cghistogram(ovi['cf',igd_ovi_comp],binsize=0.2d,min=0d,max=1d,missing=bad)
   bin3n = bin1n+bin2n
;   move 1s into top bin ...
   bin1n[4]+=bin1n[5]
   bin2n[4]+=bin2n[5]
   bin3n[4]+=bin3n[5]
  
   cgplot,[0],xran=[-0.1,1.1],yran=[0,0.7],/ysty,xticks=4,xtickv=[0.1,0.3,0.5,0.7,0.9],$
          xtickname=['0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0',' '],$
          /nodata,xminor=1,title='Covering Factor Distribution'
      
   int=0.2d
   w=0.03d
   w1=w*1.5d
   w2=w*0.5d
   z = DINDGEN(10000)*0.0001d
   c = 0.683d
   for i=0,4 do begin

      j = double(i)*int + 0.1d

      cgpolygon,[j-w1,j-w2,j-w2,j-w1,j-w1],[0,0,bin1[i],bin1[i],0],/fill,color=101
      k = bin1n[i]
      n = ctnvcomp
      Beta = IBETA(k+1,n-k+1,z)
      ill = VALUE_LOCATE(Beta,(1-c)/2)
      iul = VALUE_LOCATE(Beta,1-(1-c)/2)
      plower = z[ill]
      pupper = z[iul]
      cgoplot,[j-w],[bin1[i]],err_ylo=bin1[i]-plower,err_yhi=pupper-bin1[i]

      cgpolygon,[j-w2,j+w2,j+w2,j-w2,j-w2],[0,0,bin2[i],bin2[i],0],/fill,color=102
      k = bin2n[i]
      n = ctovicomp
      Beta = IBETA(k+1,n-k+1,z)
      ill = VALUE_LOCATE(Beta,(1-c)/2)
      iul = VALUE_LOCATE(Beta,1-(1-c)/2)
      plower = z[ill]
      pupper = z[iul]
      cgoplot,[j],[bin2[i]],err_ylo=bin2[i]-plower,err_yhi=pupper-bin2[i]

      cgpolygon,[j+w2,j+w1,j+w1,j+w2,j+w2],[0,0,bin3[i],bin3[i],0],/fill,color=103
      k = bin3n[i]
      n = ctnvcomp+ctovicomp
      Beta = IBETA(k+1,n-k+1,z)
      ill = VALUE_LOCATE(Beta,(1-c)/2)
      iul = VALUE_LOCATE(Beta,1-(1-c)/2)
      plower = z[ill]
      pupper = z[iul]
      cgoplot,[j+w],[bin3[i]],err_ylo=bin3[i]-plower,err_yhi=pupper-bin3[i]

   endfor



   al_legend,['NV','OVI','NV+OVI'],$
              color=[101,102,103],psym=[15,15,15],$
              position=[-0.1,0.7],spacing=2,/data,symsize=[2,2,2],back='White'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'hist_sig.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.5,default=4,/quiet

   bin1 = cghistogram(alog10(nv['sig',igd_nv_comp]),binsize=0.6d,min=1d,missing=bad,/freq)
   bin2 = cghistogram(alog10(ovi['sig',igd_ovi_comp]),binsize=0.6d,min=1d,missing=bad,/freq)
   bin3 = ([bin1,0]+bin2)/2d
   bin1n = cghistogram(alog10(nv['sig',igd_nv_comp]),binsize=0.6d,min=1d,missing=bad)
   bin2n = cghistogram(alog10(ovi['sig',igd_ovi_comp]),binsize=0.6d,min=1d,missing=bad)
   bin3n = [bin1n,0]+bin2n
  
   cgplot,[0],xran=[0,5],yran=[0,0.8],/ysty,xticks=3,xtickv=[1,2,3,4],$
          xtickname=['10-40','40-80','80-160','>160',' '],$
          /nodata,xminor=1,title='Log[$\sigma$ / km/s ] Distributions'
      
   int=1d
   w=0.3d
   w1=w*1.5d
   w2=w*0.5d
   z = DINDGEN(10000)*0.0001d
   c = 0.683d
   for i=0,3 do begin
      j = double(i)*int + int
      if i le n_elements(bin1)-1 then begin
         cgpolygon,[j-w1,j-w2,j-w2,j-w1,j-w1],[0,0,bin1[i],bin1[i],0],/fill,color=101
         plu = betaprob(c,bin1n[i],ctnvcomp,z)
         cgoplot,[j-w],[bin1[i]],err_ylo=bin1[i]-plu[0],err_yhi=plu[1]-bin1[i]
      endif
      if i le n_elements(bin2)-1 then begin
         cgpolygon,[j-w2,j+w2,j+w2,j-w2,j-w2],[0,0,bin2[i],bin2[i],0],/fill,color=102
         plu = betaprob(c,bin2n[i],ctovicomp,z)
         cgoplot,[j],[bin2[i]],err_ylo=bin2[i]-plu[0],err_yhi=plu[1]-bin2[i]
      endif
      if i le n_elements(bin3)-1 then begin
         cgpolygon,[j+w2,j+w1,j+w1,j+w2,j+w2],[0,0,bin3[i],bin3[i],0],/fill,color=103
         plu = betaprob(c,bin3n[i],ctnvcomp+ctovicomp,z)
         cgoplot,[j+w],[bin3[i]],err_ylo=bin3[i]-plu[0],err_yhi=plu[1]-bin3[i]
      endif
   endfor

   al_legend,['NV','OVI','NV+OVI'],$
              color=[101,102,103],psym=[15,15,15],$
              position=[0,0.8],spacing=2,/data,symsize=[2,2,2],back='White'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'hist_tau.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.5,default=4,/quiet

   bin1 = cghistogram(alog10(nv['tau',igd_nv_comp]),binsize=0.6d,min=-1.6d,missing=bad,/freq)
   bin2 = cghistogram(alog10(ovi['tau',igd_ovi_comp]),binsize=0.6d,min=-1.6d,missing=bad,/freq)
   bin3 = (bin1+bin2)/2d
   bin1n = cghistogram(alog10(nv['tau',igd_nv_comp]),binsize=0.6d,min=-1.6d,missing=bad)
   bin2n = cghistogram(alog10(ovi['tau',igd_ovi_comp]),binsize=0.6d,min=-1.6d,missing=bad)
   bin3n = bin1n+bin2n

   cgplot,[0],xran=[0,5],yran=[0,0.7],/ysty,xticks=3,xtickv=[1,2,3,4],$
      xtickname=['<0.1','0.1-0.4','0.4-1.6','1.6-5',' '],$
      /nodata,xminor=1,title='Log($\tau$$\down1243$ or $\tau$$\down1038$) Distributions'

   int=1d
   w=0.3d
   w1=w*1.5d
   w2=w*0.5d
   z = DINDGEN(10000)*0.0001d
   c = 0.683d
   for i=0,3 do begin
      j = double(i)*int + int
      if i le n_elements(bin1)-1 then begin
         cgpolygon,[j-w1,j-w2,j-w2,j-w1,j-w1],[0,0,bin1[i],bin1[i],0],/fill,color=101
         plu = betaprob(c,bin1n[i],ctnvcomp,z)
         cgoplot,[j-w],[bin1[i]],err_ylo=bin1[i]-plu[0],err_yhi=plu[1]-bin1[i]
      endif
      if i le n_elements(bin2)-1 then begin
         cgpolygon,[j-w2,j+w2,j+w2,j-w2,j-w2],[0,0,bin2[i],bin2[i],0],/fill,color=102
         plu = betaprob(c,bin2n[i],ctovicomp,z)
         cgoplot,[j],[bin2[i]],err_ylo=bin2[i]-plu[0],err_yhi=plu[1]-bin2[i]
      endif
      if i le n_elements(bin3)-1 then begin
         cgpolygon,[j+w2,j+w1,j+w1,j+w2,j+w2],[0,0,bin3[i],bin3[i],0],/fill,color=103
         plu = betaprob(c,bin3n[i],ctnvcomp+ctovicomp,z)
         cgoplot,[j+w],[bin3[i]],err_ylo=bin3[i]-plu[0],err_yhi=plu[1]-bin3[i]
      endif
   endfor

   al_legend,['NV','OVI','NV+OVI'],$
      color=[101,102,103],psym=[15,15,15],$
      position=[0,0.7],spacing=2,/data,symsize=[2,2,2],back='White'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'hist_vel.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.5,default=4,/quiet

   bin1 = cghistogram(nv['v50',igd_nv_comp],binsize=2000d,min=-8000d,missing=bad,/freq)
   bin2 = cghistogram(ovi['v50',igd_ovi_comp],binsize=2000d,min=-8000d,missing=bad,/freq)
   bin3 = (bin1+bin2)/2d
   bin1n = cghistogram(nv['v50',igd_nv_comp],binsize=2000d,min=-8000d,missing=bad)
   bin2n = cghistogram(ovi['v50',igd_ovi_comp],binsize=2000d,min=-8000d,missing=bad)
   bin3n = bin1n+bin2n

   cgplot,[0],xran=[0,6],yran=[0,0.7],/ysty,xticks=4,xtickv=[1,2,3,4,5],$
      xtickname=['-8000','-6000','-4000',$
                 '-2000','> 0',' '],$
      /nodata,xminor=1,title='Velocity (km/s) Distributions',xcharsize=1
   cgtext,0.18,0.05,'to -6000',/norm
   cgtext,0.32,0.05,'to -4000',/norm
   cgtext,0.46,0.05,'to -2000',/norm
   cgtext,0.63,0.05,'to 0',/norm
   
   int=1d
   w=0.3d
   w1=w*1.5d
   w2=w*0.5d
   z = DINDGEN(10000)*0.0001d
   c = 0.683d
   for i=0,4 do begin
      j = double(i)*int + int
      if i le n_elements(bin1)-1 then begin
         cgpolygon,[j-w1,j-w2,j-w2,j-w1,j-w1],[0,0,bin1[i],bin1[i],0],/fill,color=101
         plu = betaprob(c,bin1n[i],ctnvcomp,z)
         cgoplot,[j-w],[bin1[i]],err_ylo=bin1[i]-plu[0],err_yhi=plu[1]-bin1[i]
      endif
      if i le n_elements(bin2)-1 then begin
         cgpolygon,[j-w2,j+w2,j+w2,j-w2,j-w2],[0,0,bin2[i],bin2[i],0],/fill,color=102
         plu = betaprob(c,bin2n[i],ctovicomp,z)
         cgoplot,[j],[bin2[i]],err_ylo=bin2[i]-plu[0],err_yhi=plu[1]-bin2[i]
      endif
      if i le n_elements(bin3)-1 then begin
         cgpolygon,[j+w2,j+w1,j+w1,j+w2,j+w2],[0,0,bin3[i],bin3[i],0],/fill,color=103
         plu = betaprob(c,bin3n[i],ctnvcomp+ctovicomp,z)
         cgoplot,[j+w],[bin3[i]],err_ylo=bin3[i]-plu[0],err_yhi=plu[1]-bin3[i]
      endif
   endfor

   al_legend,['NV','OVI','NV+OVI'],$
      color=[101,102,103],psym=[15,15,15],$
      position=[0,0.7],spacing=2,/data,symsize=[2,2,2],back='White'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; TABLES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;  Detection rates

   openw,lun_tmp,tabdir+'tab_detrat.tex',/get_lun
   
   invtot = where(tabdat['nvstatus'] ne 'X',ctnvtot)
   iovitot = where(tabdat['ovistatus'] ne 'X',ctovitot)
   ianytot = where(tabdat['nvstatus'] ne 'X' $
                   OR tabdat['ovistatus'] ne 'X',ctanytot)
   ibothtot = where(tabdat['nvstatus'] ne 'X' $
                    AND tabdat['ovistatus'] ne 'X',ctbothtot)

   inv = where(tabdat['nvstatus'] eq 'F' $
               OR tabdat['nvstatus'] eq 'GF',ctnv)
   iovi = where(tabdat['ovistatus'] eq 'F' $
                OR tabdat['ovistatus'] eq 'GF',ctovi)
   iboth = where((tabdat['nvstatus'] eq 'F' $
                  OR tabdat['nvstatus'] eq 'GF') $
                  AND (tabdat['ovistatus'] eq 'F' $
                  OR tabdat['ovistatus'] eq 'GF'),ctboth)
   iany = where(tabdat['nvstatus'] eq 'F' $
                 OR tabdat['nvstatus'] eq 'GF' $
                 OR tabdat['ovistatus'] eq 'F' $
                 OR tabdat['ovistatus'] eq 'GF',ctany)

;  Total # where particular detection was possible and in particular category
   invtotlbol1 = cgsetintersection(invtot,where(tabdat['lbol'] ge 12),count=ctnvtotlbol1)
   iovitotlbol1 = cgsetintersection(iovitot,where(tabdat['lbol'] ge 12),count=ctovitotlbol1)
   ibothtotlbol1 = cgsetintersection(ibothtot,where(tabdat['lbol'] ge 12),count=ctbothtotlbol1)
   ianytotlbol1 = cgsetintersection(ianytot,where(tabdat['lbol'] ge 12),count=ctanytotlbol1)

   invtotlbol2 = cgsetintersection(invtot,where(tabdat['lbol'] lt 12),count=ctnvtotlbol2)
   iovitotlbol2 = cgsetintersection(iovitot,where(tabdat['lbol'] lt 12),count=ctovitotlbol2)
   ibothtotlbol2 = cgsetintersection(ibothtot,where(tabdat['lbol'] lt 12),count=ctbothtotlbol2)
   ianytotlbol2 = cgsetintersection(ianytot,where(tabdat['lbol'] lt 12),count=ctanytotlbol2)

;  this assumes all X-ray measurements for a source are either N_H = 0 OR N_H > 0
   invtotnhnz = cgsetintersection(invtot,where((tabdat['nhxray',*,0] gt 22d AND $
      tabdat['nhxray',*,0] ne bad) OR (tabdat['nhxray_lim',*,0] gt 22d AND $
      tabdat['nhxray_lim',*,0] ne bad)),count=ctnvtotnhnz)
   iovitotnhnz = cgsetintersection(iovitot,where((tabdat['nhxray',*,0] gt 22d AND $
      tabdat['nhxray',*,0] ne bad) OR (tabdat['nhxray_lim',*,0] gt 22d AND $
      tabdat['nhxray_lim',*,0] ne bad)),count=ctovitotnhnz)
   ibothtotnhnz = cgsetintersection(ibothtot,where((tabdat['nhxray',*,0] gt 22d AND $
      tabdat['nhxray',*,0] ne bad) OR (tabdat['nhxray_lim',*,0] gt 22d AND $
      tabdat['nhxray_lim',*,0] ne bad)),count=ctbothtotnhnz)
   ianytotnhnz = cgsetintersection(ianytot,where((tabdat['nhxray',*,0] gt 22d AND $
      tabdat['nhxray',*,0] ne bad) OR (tabdat['nhxray_lim',*,0] gt 22d AND $
      tabdat['nhxray_lim',*,0] ne bad)),count=ctanytotnhnz)

   invtotnhz = cgsetintersection(invtot,where(tabdat['nhxray',*,0] le 22d OR $
      tabdat['nhxray_lim',*,0] le 22d),count=ctnvtotnhz)
   iovitotnhz = cgsetintersection(iovitot,where(tabdat['nhxray',*,0] le 22d OR $
      tabdat['nhxray_lim',*,0] le 22d),count=ctovitotnhz)
   ibothtotnhz = cgsetintersection(ibothtot,where(tabdat['nhxray',*,0] le 22d OR $
      tabdat['nhxray_lim',*,0] le 22d),count=ctbothtotnhz)
   ianytotnhz = cgsetintersection(ianytot,where(tabdat['nhxray',*,0] le 22d OR $
      tabdat['nhxray_lim',*,0] le 22d),count=ctanytotnhz)

   invtotaox1 = cgsetintersection(invtot,where(tabdat['alphaox'] ge -1.6$
      AND tabdat['alphaox'] ne bad),count=ctnvtotaox1)
   iovitotaox1 = cgsetintersection(iovitot,where(tabdat['alphaox'] ge -1.6$
      AND tabdat['alphaox'] ne bad),count=ctovitotaox1)
   ibothtotaox1 = cgsetintersection(ibothtot,where(tabdat['alphaox'] ge -1.6$
      AND tabdat['alphaox'] ne bad),count=ctbothtotaox1)
   ianytotaox1 = cgsetintersection(ianytot,where(tabdat['alphaox'] ge -1.6$
      AND tabdat['alphaox'] ne bad),count=ctanytotaox1)

   invtotaox2 = cgsetintersection(invtot,where(tabdat['alphaox'] lt -1.6),count=ctnvtotaox2)
   iovitotaox2 = cgsetintersection(iovitot,where(tabdat['alphaox'] lt -1.6),count=ctovitotaox2)
   ibothtotaox2 = cgsetintersection(ibothtot,where(tabdat['alphaox'] lt -1.6),count=ctbothtotaox2)
   ianytotaox2 = cgsetintersection(ianytot,where(tabdat['alphaox'] lt -1.6),count=ctanytotaox2)

   invtotagnfrac1 = cgsetintersection(invtot,where(tabdat['agnfrac'] ge 0.95),count=ctnvtotagnfrac1)
   iovitotagnfrac1 = cgsetintersection(iovitot,where(tabdat['agnfrac'] ge 0.95),count=ctovitotagnfrac1)
   ibothtotagnfrac1 = cgsetintersection(ibothtot,where(tabdat['agnfrac'] ge 0.95),count=ctbothtotagnfrac1)
   ianytotagnfrac1 = cgsetintersection(ianytot,where(tabdat['agnfrac'] ge 0.95),count=ctanytotagnfrac1)

   invtotagnfrac2 = cgsetintersection(invtot,where(tabdat['agnfrac'] lt 0.95),count=ctnvtotagnfrac2)
   iovitotagnfrac2 = cgsetintersection(iovitot,where(tabdat['agnfrac'] lt 0.95),count=ctovitotagnfrac2)
   ibothtotagnfrac2 = cgsetintersection(ibothtot,where(tabdat['agnfrac'] lt 0.95),count=ctbothtotagnfrac2)
   ianytotagnfrac2 = cgsetintersection(ianytot,where(tabdat['agnfrac'] lt 0.95),count=ctanytotagnfrac2)

   invtotluv1 = cgsetintersection(invtot,where(tabdat['luv'] ge 45.1),count=ctnvtotluv1)
   iovitotluv1 = cgsetintersection(iovitot,where(tabdat['luv'] ge 45.1),count=ctovitotluv1)
   ibothtotluv1 = cgsetintersection(ibothtot,where(tabdat['luv'] ge 45.1),count=ctbothtotluv1)
   ianytotluv1 = cgsetintersection(ianytot,where(tabdat['luv'] ge 45.1),count=ctanytotluv1)

   invtotluv2 = cgsetintersection(invtot,where(tabdat['luv'] lt 45.1),count=ctnvtotluv2)
   iovitotluv2 = cgsetintersection(iovitot,where(tabdat['luv'] lt 45.1),count=ctovitotluv2)
   ibothtotluv2 = cgsetintersection(ibothtot,where(tabdat['luv'] lt 45.1),count=ctbothtotluv2)
   ianytotluv2 = cgsetintersection(ianytot,where(tabdat['luv'] lt 45.1),count=ctanytotluv2)

;  Total # where particular detection was made and in particular category
   invlbol1 = cgsetintersection(inv,where(tabdat['lbol'] ge 12),count=ctnvlbol1)
   iovilbol1 = cgsetintersection(iovi,where(tabdat['lbol'] ge 12),count=ctovilbol1)
   ibothlbol1 = cgsetintersection(iboth,where(tabdat['lbol'] ge 12),count=ctbothlbol1)
   ianylbol1 = cgsetintersection(iany,where(tabdat['lbol'] ge 12),count=ctanylbol1)

   invlbol2 = cgsetintersection(inv,where(tabdat['lbol'] lt 12),count=ctnvlbol2)
   iovilbol2 = cgsetintersection(iovi,where(tabdat['lbol'] lt 12),count=ctovilbol2)
   ibothlbol2 = cgsetintersection(iboth,where(tabdat['lbol'] lt 12),count=ctbothlbol2)
   ianylbol2 = cgsetintersection(iany,where(tabdat['lbol'] lt 12),count=ctanylbol2)

   invnhnz = cgsetintersection(inv,where((tabdat['nhxray',*,0] gt 22d AND $
      tabdat['nhxray',*,0] ne bad) OR (tabdat['nhxray_lim',*,0] gt 22d AND $
      tabdat['nhxray_lim',*,0] ne bad)),count=ctnvnhnz)
   iovinhnz = cgsetintersection(iovi,where((tabdat['nhxray',*,0] gt 22d AND $
      tabdat['nhxray',*,0] ne bad) OR (tabdat['nhxray_lim',*,0] gt 22d AND $
      tabdat['nhxray_lim',*,0] ne bad)),count=ctovinhnz)
   ibothnhnz = cgsetintersection(iboth,where((tabdat['nhxray',*,0] gt 22d AND $
      tabdat['nhxray',*,0] ne bad) OR (tabdat['nhxray_lim',*,0] gt 22d AND $
      tabdat['nhxray_lim',*,0] ne bad)),count=ctbothnhnz)
   ianynhnz = cgsetintersection(iany,where((tabdat['nhxray',*,0] gt 22d AND $
      tabdat['nhxray',*,0] ne bad) OR (tabdat['nhxray_lim',*,0] gt 22d AND $
      tabdat['nhxray_lim',*,0] ne bad)),count=ctanynhnz)

   invnhz = cgsetintersection(inv,where(tabdat['nhxray',*,0] le 22d OR $
      tabdat['nhxray_lim',*,0] le 22d),count=ctnvnhz)
   iovinhz = cgsetintersection(iovi,where(tabdat['nhxray',*,0] le 22d OR $
      tabdat['nhxray_lim',*,0] le 22d),count=ctovinhz)
   ibothnhz = cgsetintersection(iboth,where(tabdat['nhxray',*,0] le 22d OR $
      tabdat['nhxray_lim',*,0] le 22d),count=ctbothnhz)
   ianynhz = cgsetintersection(iany,where(tabdat['nhxray',*,0] le 22d OR $
      tabdat['nhxray_lim',*,0] le 22d),count=ctanynhz)

   invaox1 = cgsetintersection(inv,where(tabdat['alphaox'] ge -1.6$
      AND tabdat['alphaox'] ne bad),count=ctnvaox1)
   ioviaox1 = cgsetintersection(iovi,where(tabdat['alphaox'] ge -1.6$
      AND tabdat['alphaox'] ne bad),count=ctoviaox1)
   ibothaox1 = cgsetintersection(iboth,where(tabdat['alphaox'] ge -1.6$
      AND tabdat['alphaox'] ne bad),count=ctbothaox1)
   ianyaox1 = cgsetintersection(iany,where(tabdat['alphaox'] ge -1.6$
      AND tabdat['alphaox'] ne bad),count=ctanyaox1)

   invaox2 = cgsetintersection(inv,where(tabdat['alphaox'] lt -1.6),count=ctnvaox2)
   ioviaox2 = cgsetintersection(iovi,where(tabdat['alphaox'] lt -1.6),count=ctoviaox2)
   ibothaox2 = cgsetintersection(iboth,where(tabdat['alphaox'] lt -1.6),count=ctbothaox2)
   ianyaox2 = cgsetintersection(iany,where(tabdat['alphaox'] lt -1.6),count=ctanyaox2)

   invagnfrac1 = cgsetintersection(inv,where(tabdat['agnfrac'] ge 0.95),count=ctnvagnfrac1)
   ioviagnfrac1 = cgsetintersection(iovi,where(tabdat['agnfrac'] ge 0.95),count=ctoviagnfrac1)
   ibothagnfrac1 = cgsetintersection(iboth,where(tabdat['agnfrac'] ge 0.95),count=ctbothagnfrac1)
   ianyagnfrac1 = cgsetintersection(iany,where(tabdat['agnfrac'] ge 0.95),count=ctanyagnfrac1)

   invagnfrac2 = cgsetintersection(inv,where(tabdat['agnfrac'] lt 0.95),count=ctnvagnfrac2)
   ioviagnfrac2 = cgsetintersection(iovi,where(tabdat['agnfrac'] lt 0.95),count=ctoviagnfrac2)
   ibothagnfrac2 = cgsetintersection(iboth,where(tabdat['agnfrac'] lt 0.95),count=ctbothagnfrac2)
   ianyagnfrac2 = cgsetintersection(iany,where(tabdat['agnfrac'] lt 0.95),count=ctanyagnfrac2)

   invluv1 = cgsetintersection(inv,where(tabdat['luv'] ge 45.1),count=ctnvluv1)
   ioviluv1 = cgsetintersection(iovi,where(tabdat['luv'] ge 45.1),count=ctoviluv1)
   ibothluv1 = cgsetintersection(iboth,where(tabdat['luv'] ge 45.1),count=ctbothluv1)
   ianyluv1 = cgsetintersection(iany,where(tabdat['luv'] ge 45.1),count=ctanyluv1)

   invluv2 = cgsetintersection(inv,where(tabdat['luv'] lt 45.1),count=ctnvluv2)
   ioviluv2 = cgsetintersection(iovi,where(tabdat['luv'] lt 45.1),count=ctoviluv2)
   ibothluv2 = cgsetintersection(iboth,where(tabdat['luv'] lt 45.1),count=ctbothluv2)
   ianyluv2 = cgsetintersection(iany,where(tabdat['luv'] lt 45.1),count=ctanyluv2)

   z = DINDGEN(10000)*0.0001d
   c = 0.683d
   ;c = 0.954d

;  All Quasars

   printf,lun_tmp,'\cutinhead{All Quasars}'
   plu = betaprob(c,ctnv,ctnvtot,z)
   printf,lun_tmp,'N V',amp,ctnv,amp,ctnvtot,amp,double(ctnv)/double(ctnvtot),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctovi,ctovitot,z)
   printf,lun_tmp,'O VI',amp,ctovi,amp,ctovitot,amp,double(ctovi)/double(ctovitot),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctboth,ctbothtot,z)
   printf,lun_tmp,'Both',amp,ctboth,amp,ctbothtot,amp,double(ctboth)/double(ctbothtot),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctany,ctanytot,z)
   printf,lun_tmp,'Any',amp,ctany,amp,ctanytot,amp,double(ctany)/double(ctanytot),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'


;  Lbol >= 12.0

   printf,lun_tmp,'\cutinhead{log $L_{\rm BOL}/L_\odot$ $\geq$ 12.0}'
   plu = betaprob(c,ctnvlbol1,ctnvtotlbol1,z)
   printf,lun_tmp,'N V',amp,ctnvlbol1,amp,ctnvtotlbol1,amp,$
          double(ctnvlbol1)/double(ctnvtotlbol1),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctovilbol1,ctovitotlbol1,z)
   printf,lun_tmp,'O VI',amp,ctovilbol1,amp,ctovitotlbol1,amp,$
          double(ctovilbol1)/double(ctovitotlbol1),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctbothlbol1,ctbothtotlbol1,z)
   printf,lun_tmp,'Both',amp,ctbothlbol1,amp,ctbothtotlbol1,amp,$
          double(ctbothlbol1)/double(ctbothtotlbol1),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctanylbol1,ctanytotlbol1,z)
   printf,lun_tmp,'Any',amp,ctanylbol1,amp,ctanytotlbol1,amp,$
          double(ctanylbol1)/double(ctanytotlbol1),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'

;  Lbol < 12.0

   printf,lun_tmp,'\cutinhead{log $L_{\rm BOL}/L_\odot$ $<$ 12.0}'
   plu = betaprob(c,ctnvlbol2,ctnvtotlbol2,z)
   printf,lun_tmp,'N V',amp,ctnvlbol2,amp,ctnvtotlbol2,amp,$
          double(ctnvlbol2)/double(ctnvtotlbol2),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctovilbol2,ctovitotlbol2,z)
   printf,lun_tmp,'O VI',amp,ctovilbol2,amp,ctovitotlbol2,amp,$
          double(ctovilbol2)/double(ctovitotlbol2),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctbothlbol2,ctbothtotlbol2,z)
   printf,lun_tmp,'Both',amp,ctbothlbol2,amp,ctbothtotlbol2,amp,$
          double(ctbothlbol2)/double(ctbothtotlbol2),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctanylbol2,ctanytotlbol2,z)
   printf,lun_tmp,'Any',amp,ctanylbol2,amp,ctanytotlbol2,amp,$
          double(ctanylbol2)/double(ctanytotlbol2),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'


;;  L1125 >= 45.1
;
;   printf,lun_tmp,'\cutinhead{log $L_{\rm 1125} $\geq$ 45.1}'
;   plu = betaprob(c,ctnvluv1,ctnvtotluv1,z)
;   printf,lun_tmp,'N V',amp,ctnvluv1,amp,ctnvtotluv1,amp,$
;          double(ctnvluv1)/double(ctnvtotluv1),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctoviluv1,ctovitotluv1,z)
;   printf,lun_tmp,'O VI',amp,ctoviluv1,amp,ctovitotluv1,amp,$
;          double(ctoviluv1)/double(ctovitotluv1),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctbothluv1,ctbothtotluv1,z)
;   printf,lun_tmp,'Both',amp,ctbothluv1,amp,ctbothtotluv1,amp,$
;          double(ctbothluv1)/double(ctbothtotluv1),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctanyluv1,ctanytotluv1,z)
;   printf,lun_tmp,'Any',amp,ctanyluv1,amp,ctanytotluv1,amp,$
;          double(ctanyluv1)/double(ctanytotluv1),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;
;;  L1125 < 45.1
;
;   printf,lun_tmp,'\cutinhead{log $L_{\rm 1125}$ $<$ 45.1}'
;   plu = betaprob(c,ctnvluv2,ctnvtotluv2,z)
;   printf,lun_tmp,'N V',amp,ctnvluv2,amp,ctnvtotluv2,amp,$
;          double(ctnvluv2)/double(ctnvtotluv2),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctoviluv2,ctovitotluv2,z)
;   printf,lun_tmp,'O VI',amp,ctoviluv2,amp,ctovitotluv2,amp,$
;          double(ctoviluv2)/double(ctovitotluv2),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctbothluv2,ctbothtotluv2,z)
;   printf,lun_tmp,'Both',amp,ctbothluv2,amp,ctbothtotluv2,amp,$
;          double(ctbothluv2)/double(ctbothtotluv2),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctanyluv2,ctanytotluv2,z)
;   printf,lun_tmp,'Any',amp,ctanyluv2,amp,ctanytotluv2,amp,$
;          double(ctanyluv2)/double(ctanytotluv2),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'

;  N_H > 0

   printf,lun_tmp,'\cutinhead{$N_{\rm H}$ $>$ $10^{22}$ cm$^{-2}$}'
   plu = betaprob(c,ctnvnhnz,ctnvtotnhnz,z)
   printf,lun_tmp,'N V',amp,ctnvnhnz,amp,ctnvtotnhnz,amp,$
          double(ctnvnhnz)/double(ctnvtotnhnz),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctovinhnz,ctovitotnhnz,z)
   printf,lun_tmp,'O VI',amp,ctovinhnz,amp,ctovitotnhnz,amp,$
          double(ctovinhnz)/double(ctovitotnhnz),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctbothnhnz,ctbothtotnhnz,z)
   printf,lun_tmp,'Both',amp,ctbothnhnz,amp,ctbothtotnhnz,amp,$
          double(ctbothnhnz)/double(ctbothtotnhnz),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctanynhnz,ctanytotnhnz,z)
   printf,lun_tmp,'Any',amp,ctanynhnz,amp,ctanytotnhnz,amp,$
          double(ctanynhnz)/double(ctanytotnhnz),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'

;  N_H = 0

   printf,lun_tmp,'\cutinhead{$N_{\rm H}$ $\leq$ $10^{22}$ cm$^{-2}$}'
   plu = betaprob(c,ctnvnhz,ctnvtotnhz,z)
   printf,lun_tmp,'N V',amp,ctnvnhz,amp,ctnvtotnhz,amp,$
          double(ctnvnhz)/double(ctnvtotnhz),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctovinhz,ctovitotnhz,z)
   printf,lun_tmp,'O VI',amp,ctovinhz,amp,ctovitotnhz,amp,$
          double(ctovinhz)/double(ctovitotnhz),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctbothnhz,ctbothtotnhz,z)
   printf,lun_tmp,'Both',amp,ctbothnhz,amp,ctbothtotnhz,amp,$
          double(ctbothnhz)/double(ctbothtotnhz),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctanynhz,ctanytotnhz,z)
   printf,lun_tmp,'Any',amp,ctanynhz,amp,ctanytotnhz,amp,$
          double(ctanynhz)/double(ctanytotnhz),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'

;  alphaox >= -1.6

   printf,lun_tmp,'\cutinhead{$\alpha_{ox}$ $\geq$ $-$1.6}'
   plu = betaprob(c,ctnvaox1,ctnvtotaox1,z)
   printf,lun_tmp,'N V',amp,ctnvaox1,amp,ctnvtotaox1,amp,$
          double(ctnvaox1)/double(ctnvtotaox1),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctoviaox1,ctovitotaox1,z)
   printf,lun_tmp,'O VI',amp,ctoviaox1,amp,ctovitotaox1,amp,$
          double(ctoviaox1)/double(ctovitotaox1),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctbothaox1,ctbothtotaox1,z)
   printf,lun_tmp,'Both',amp,ctbothaox1,amp,ctbothtotaox1,amp,$
          double(ctbothaox1)/double(ctbothtotaox1),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctanyaox1,ctanytotaox1,z)
   printf,lun_tmp,'Any',amp,ctanyaox1,amp,ctanytotaox1,amp,$
          double(ctanyaox1)/double(ctanytotaox1),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'

;  alphaox < -1.6

   printf,lun_tmp,'\cutinhead{$\alpha_{ox}$ $<$ $-$1.6}'
   plu = betaprob(c,ctnvaox2,ctnvtotaox2,z)
   printf,lun_tmp,'N V',amp,ctnvaox2,amp,ctnvtotaox2,amp,$
          double(ctnvaox2)/double(ctnvtotaox2),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctoviaox2,ctovitotaox2,z)
   printf,lun_tmp,'O VI',amp,ctoviaox2,amp,ctovitotaox2,amp,$
          double(ctoviaox2)/double(ctovitotaox2),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctbothaox2,ctbothtotaox2,z)
   printf,lun_tmp,'Both',amp,ctbothaox2,amp,ctbothtotaox2,amp,$
          double(ctbothaox2)/double(ctbothtotaox2),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
   plu = betaprob(c,ctanyaox2,ctanytotaox2,z)
   printf,lun_tmp,'Any',amp,ctanyaox2,amp,ctanytotaox2,amp,$
          double(ctanyaox2)/double(ctanytotaox2),$
          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'


;;  agnfrac >= 0.95
;
;   printf,lun_tmp,'\cutinhead{AGN fraction $\geq$ 0.95}'
;   plu = betaprob(c,ctnvagnfrac1,ctnvtotagnfrac1,z)
;   printf,lun_tmp,'N V',amp,ctnvagnfrac1,amp,ctnvtotagnfrac1,amp,$
;          double(ctnvagnfrac1)/double(ctnvtotagnfrac1),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctoviagnfrac1,ctovitotagnfrac1,z)
;   printf,lun_tmp,'O VI',amp,ctoviagnfrac1,amp,ctovitotagnfrac1,amp,$
;          double(ctoviagnfrac1)/double(ctovitotagnfrac1),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctbothagnfrac1,ctbothtotagnfrac1,z)
;   printf,lun_tmp,'Both',amp,ctbothagnfrac1,amp,ctbothtotagnfrac1,amp,$
;          double(ctbothagnfrac1)/double(ctbothtotagnfrac1),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctanyagnfrac1,ctanytotagnfrac1,z)
;   printf,lun_tmp,'Any',amp,ctanyagnfrac1,amp,ctanytotagnfrac1,amp,$
;          double(ctanyagnfrac1)/double(ctanytotagnfrac1),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;
;;  agnfrac < 0.95
;
;   printf,lun_tmp,'\cutinhead{AGN fraction $<$ 0.95}'
;   plu = betaprob(c,ctnvagnfrac2,ctnvtotagnfrac2,z)
;   printf,lun_tmp,'N V',amp,ctnvagnfrac2,amp,ctnvtotagnfrac2,amp,$
;          double(ctnvagnfrac2)/double(ctnvtotagnfrac2),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctoviagnfrac2,ctovitotagnfrac2,z)
;   printf,lun_tmp,'O VI',amp,ctoviagnfrac2,amp,ctovitotagnfrac2,amp,$
;          double(ctoviagnfrac2)/double(ctovitotagnfrac2),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctbothagnfrac2,ctbothtotagnfrac2,z)
;   printf,lun_tmp,'Both',amp,ctbothagnfrac2,amp,ctbothtotagnfrac2,amp,$
;          double(ctbothagnfrac2)/double(ctbothtotagnfrac2),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'
;   plu = betaprob(c,ctanyagnfrac2,ctanytotagnfrac2,z)
;   printf,lun_tmp,'Any',amp,ctanyagnfrac2,amp,ctanytotagnfrac2,amp,$
;          double(ctanyagnfrac2)/double(ctanytotagnfrac2),$
;          ' (',plu[0],' $-$ ',plu[1],')',dslash,$
;          format='(A5,A3,I3,A3,I3,A3,D5.2,A0,D4.2,A0,D4.2,A0,A0)'

   free_lun,lun_tmp

; Fit results

   openw,lun_tmp,tabdir+'tab_fitresults.tex',/get_lun
   
   for i=0,ncos-1 do begin
      
      labprint=0b
      if tabdat['nvstatus',i] eq 'F' OR tabdat['nvstatus',i] eq 'GF' then begin
         printf,lun_tmp,tabdat['gal',i],amp,'N~V',amp,nv['weq_A',i,0],$
            lerr0,nv['weq_A',i,1],lerr1,nv['weq_A',i,2],lerr2,amp,$
            nv['vwtavg',i,0],lerr0,nv['vwtavg',i,1],lerr1,nv['vwtavg',i,2],lerr2,$
            amp,nv['vwtrms',i],lerr0,nv['vwtrms',i,1],lerr1,nv['vwtrms',i,2],lerr2,$
            amp,nv['ncomp',i],dslash,$
            format='(A12,A3,A5,A3,'+$
            'D6.2,A0,D-5.3,A0,D-5.3,A0,A3,'+$
            'D7.1,A0,D-5.1,A0,D-5.1,A0,A3,'+$
            'D7.1,A0,D-5.1,A0,D-5.1,A0,A3,'+$
            'I3,A3)'
         labprint=1b
      endif else if nv['weq_lim_A',i] ne bad then begin
         if labprint then labtmp = '' else labtmp = tabdat['gal',i]
         printf,lun_tmp,labtmp,amp,'N~V',amp,'$<$',nv['weq_lim_A',i],$
            amp,amp,amp,dslash,$
            format='(A12,A3,A5,A3,A3,D0.2,4A3)'
         labprint=1b
      endif
      if tabdat['ovistatus',i] eq 'F' OR tabdat['ovistatus',i] eq 'GF' then begin
         if labprint then labtmp = '' else labtmp = tabdat['gal',i]
         printf,lun_tmp,labtmp,amp,'O~VI',amp,ovi['weq_A',i,0],$
            lerr0,ovi['weq_A',i,1],lerr1,ovi['weq_A',i,2],lerr2,amp,$
            ovi['vwtavg',i,0],lerr0,ovi['vwtavg',i,1],lerr1,ovi['vwtavg',i,2],lerr2,$
            amp,ovi['vwtrms',i],lerr0,ovi['vwtrms',i,1],lerr1,ovi['vwtrms',i,2],lerr2,$
            amp,ovi['ncomp',i],dslash,$
            format='(A12,A3,A5,A3,'+$
            'D6.2,A0,D-5.3,A0,D-5.3,A0,A3,'+$
            'D7.1,A0,D-5.1,A0,D-5.1,A0,A3,'+$
            'D7.1,A0,D-5.1,A0,D-5.1,A0,A3,'+$
            'I3,A3)'
         labprint=1b
      endif else if ovi['weq_lim_A',i] ne bad then begin
         if labprint then labtmp = '' else labtmp = tabdat['gal',i]
         printf,lun_tmp,labtmp,amp,'O~VI',amp,'$<$',ovi['weq_lim_A',i],$
            amp,amp,amp,dslash,$
            format='(A12,A3,A5,A3,A3,D0.2,4A3)'
         labprint=1b
      endif

      if pv['ncomp',i] ne 0 then begin
         if labprint then labtmp = '' else labtmp = tabdat['gal',i]
         printf,lun_tmp,labtmp,amp,'P~V',amp,pv['weq_A',i,0],$
            lerr0,pv['weq_A',i,1],lerr1,pv['weq_A',i,2],lerr2,amp,$
            pv['vwtavg',i,0],lerr0,pv['vwtavg',i,1],lerr1,pv['vwtavg',i,2],lerr2,$
            amp,pv['vwtrms',i],lerr0,pv['vwtrms',i,1],lerr1,pv['vwtrms',i,2],lerr2,$
            amp,pv['ncomp',i],dslash,$
            format='(A12,A3,A5,A3,'+$
            'D6.2,A0,D-5.3,A0,D-5.3,A0,A3,'+$
            'D7.1,A0,D-5.1,A0,D-5.1,A0,A3,'+$
            'D7.1,A0,D-5.1,A0,D-5.1,A0,A3,'+$
            'I3,A3)'
      endif
      
   endfor
   
   free_lun,lun_tmp

   ; Fit results

   openw,lun_tmp,tabdir+'tab_table1quants.tex',/get_lun

   printf,lun_tmp,'#Col 2: log(Lbol/Lsun) [corrected to cosmology for paper]
   printf,lun_tmp,'#Col 3: log(LIR/Lsun) [corrected to cosmology for paper]
   printf,lun_tmp,'#Col 4: log(lam L_lam) at 1125 A; 1.5% bandpass; corrected for Gal. extinction
   printf,lun_tmp,'#Col 5: AGNfraction_error_low^error_hi
   printf,lun_tmp,'#Col 6: log(MBH/Msun)_error_low^error_hi
   printf,lun_tmp,'#Col 7: log(Edd.rat.)_error_low^error_hi
   for i=0,ncos-1 do begin
      if tabdatagnfrac[i] eq bad then $
         agnfracstring = '1\tablenotemark{a}' $
      else $
         agnfracstring = string(tabdat['agnfrac',i],'$_{-',$
         tabdat['agnfrac',i]-tabdat['agnfraclb',i],'}^{+',$
         tabdat['agnfracub',i]-tabdat['agnfrac',i],'}$',$
         format='(D0.3,A0,D0.3,A0,D0.3,A0)')
      printf,lun_tmp,tabdat['gal',i],tabdat['lbol',i],tabdat['lir',i],$
         tabdat['luv',i],$
         agnfracstring,$
         string(lmbh[i],'$_{-',lmbh_errlo[i],'}^{+',$
         lmbh_errhi[i],'}$',format='(D0.2,A0,D0.2,A0,D0.2,A0)'),$
         string(leddrat[i],'$_{-',leddrat_errlo[i],'}^{+',$
         leddrat_errhi[i],'}$',format='(D0.2,A0,D0.2,A0,D0.2,A0)'),$
         format='(A-12,3D8.2,A28,2A25)'
   endfor

   free_lun,lun_tmp

;  free up stats tables
   free_lun,lun_stat
   free_lun,lun_reg_tex

end
