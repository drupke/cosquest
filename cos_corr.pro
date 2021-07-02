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

   qsotab='/Users/drupke/Box Sync/qsos/qsos.csv'
   cosdir='/Users/drupke/Box Sync/cosquest/'
   fitdir=cosdir+'fits/'
   plotdir=cosdir+'plots/correlations/'
   tabdir=cosdir+'tables/'

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

   
;  X-ray data requires special handling because of multiple measurements
   tabdat['n_xray'] = dblarr(ncos)
   tabdat['gamxray'] = dblarr(ncos,6)+bad
   tabdat['gamxray_errlo'] = dblarr(ncos,6)+bad
   tabdat['gamxray_errhi'] = dblarr(ncos,6)+bad
   tabdat['nhxray'] = dblarr(ncos,6)+bad
   tabdat['nhxray_errlo'] = dblarr(ncos,6)+bad
   tabdat['nhxray_errhi'] = dblarr(ncos,6)+bad
   tabdat['fsoftxray'] = dblarr(ncos,6)+bad
   tabdat['fsoftxray_errlo'] = dblarr(ncos,6)+bad
   tabdat['fsoftxray_errhi'] = dblarr(ncos,6)+bad
   tabdat['fhardxray'] = dblarr(ncos,6)+bad
   tabdat['fhardxray_errlo'] = dblarr(ncos,6)+bad
   tabdat['fhardxray_errhi'] = dblarr(ncos,6)+bad
   tabdat['lsoftxray'] = dblarr(ncos,6)+bad
   tabdat['lhardxray'] = dblarr(ncos,6)+bad
   tabdat['ltotxray'] = dblarr(ncos,6)+bad
   tabdat['lsoftratxray'] = dblarr(ncos,6)+bad
   tabdat['lxlbol'] = dblarr(ncos,6)+bad
   j = -1 ; galaxy index
   k = 0 ; zero x-ray component index
   for i=0,nlines-1 do begin
      if gal[i] ne '' AND cossamp[i] eq 'V18' then begin
         j++ ; increment galaxy index
         if gamxray[i] ne bad then begin
            k = 0 ; re-zero x-ray component index
;           record first x-ray measurement
            tabdat['n_xray',j] = k+1
            tabdat['gamxray',j,k] = gamxray[i]
            tabdat['gamxray_errlo',j,k] = gamxray_errlo[i]
            tabdat['gamxray_errhi',j,k] = gamxray_errhi[i]
            tabdat['nhxray',j,k] = nhxray[i]
            tabdat['nhxray_errlo',j,k] = nhxray_errlo[i]
            tabdat['nhxray_errhi',j,k] = nhxray_errhi[i]
            tabdat['fhardxray',j,k] = fhardxray[i]
            tabdat['fhardxray_errlo',j,k] = fhardxray_errlo[i]
            tabdat['fhardxray_errhi',j,k] = fhardxray_errhi[i]
            if fsoftxray[i] ne bad then begin
               tabdat['lsoftxray',j,k] = lsoftxray[i]
               tabdat['lhardxray',j,k] = lhardxray[i]
               tabdat['ltotxray',j,k] = lsoftxray[i] + lhardxray[i]
               tabdat['lsoftratxray',j,k] = lsoftxray[i] / tabdat['ltotxray',j,k]
               tabdat['lxlbol',j,k] = $
                  alog10(tabdat['ltotxray',j,k])+44d - llsun_ergps - tabdat['lbol',j]
               tabdat['fsoftxray',j,k] = fsoftxray[i]
               tabdat['fsoftxray_errlo',j,k] = fsoftxray_errlo[i]
               tabdat['fsoftxray_errhi',j,k] = fsoftxray_errhi[i]
            endif
            k++
         endif else begin
            tabdat['n_xray',j] = 0 ; record number of x-ray components
            k = 0 ; re-zero x-ray component index
         endelse
      endif else begin
         if k gt 0 then begin
            if gamxray[i] ne bad then begin
               tabdat['n_xray',j] = k+1
               tabdat['gamxray',j,k] = gamxray[i]
               tabdat['gamxray_errlo',j,k] = gamxray_errlo[i]
               tabdat['gamxray_errhi',j,k] = gamxray_errhi[i]
               tabdat['nhxray',j,k] = nhxray[i]
               tabdat['nhxray_errlo',j,k] = nhxray_errlo[i]
               tabdat['nhxray_errhi',j,k] = nhxray_errhi[i]
               tabdat['fhardxray',j,k] = fhardxray[i]
               tabdat['fhardxray_errlo',j,k] = fhardxray_errlo[i]
               tabdat['fhardxray_errhi',j,k] = fhardxray_errhi[i]
               tabdat['lhardxray',j,k] = lhardxray[i]
               if fsoftxray[i] ne bad then begin
                  tabdat['lsoftxray',j,k] = lsoftxray[i]
                  tabdat['fsoftxray',j,k] = fsoftxray[i]
                  tabdat['fsoftxray_errlo',j,k] = fsoftxray_errlo[i]
                  tabdat['fsoftxray_errhi',j,k] = fsoftxray_errhi[i]
                  tabdat['ltotxray',j,k] = lsoftxray[i] + lhardxray[i]
                  tabdat['lsoftratxray',j,k] = lsoftxray[i] / tabdat['ltotxray',j,k]
                  tabdat['lxlbol',j,k] = $
                     alog10(tabdat['ltotxray',j,k])+44d - llsun_ergps - tabdat['lbol',j]
               endif
               k++
            endif else begin
               k = 0 ; re-zero x-ray component index
            endelse
         endif
      endelse
   endfor

; Compute physical quantities

;  Median alpha_ox
   igdalphaox = where(tabdat['alphaox'] ne bad)
   print,'Median alpha_ox: ',median(tabdat['alphaox',igdalphaox]),format='(A0,D0.2)'




;  Set AGN fraction to 1 if we don't have a measurement
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
   igdlbol = where(tabdat['lbol'] ne bad)
   lagn[igdlbol] = tabdat['lbol',igdlbol] + $
                   alog10(tabdat['agnfrac',igdlbol])
   lagnlb[igdlbol] = tabdat['lbol',igdlbol] + $
                     alog10(tabdat['agnfraclb',igdlbol])
   lagnub[igdlbol] = tabdat['lbol',igdlbol] + $
                     alog10(tabdat['agnfracub',igdlbol])
                     
;  Black hole masses
   lmbh = dblarr(ncos)+bad
   lmbh_errlo = dblarr(ncos)+bad
   lmbh_errhi = dblarr(ncos)+bad
   igdrev = where(tabdat['lmbh_rev'] ne bad)
   igdrevrm = where(tabdat['lmbh_rev'] ne bad AND $
                    tabdat['lmbh_rev_type'] eq 'RM')
   igdrevse = where(tabdat['lmbh_rev'] ne bad AND $
                    tabdat['lmbh_rev_type'] eq 'SE')
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
   nv['weq'] = dblarr(ncos)+bad
   nv['weq_A'] = dblarr(ncos)+bad
   nv['vwtavg'] = dblarr(ncos)+bad
   nv['vwtrms'] = dblarr(ncos)+bad
   nv['v50'] = dblarr(ncos,maxncomp)+bad
   nv['cf'] = dblarr(ncos,maxncomp)+bad
   nv['tau'] = dblarr(ncos,maxncomp)+bad
   nv['sig'] = dblarr(ncos,maxncomp)+bad
   ovi = orderedhash()
   ovi['ncomp'] = intarr(ncos)
   ovi['weq'] = dblarr(ncos)+bad
   ovi['weq_A'] = dblarr(ncos)+bad
   ovi['vwtavg'] = dblarr(ncos)+bad
   ovi['vwtrms'] = dblarr(ncos)+bad
   ovi['v50'] = dblarr(ncos,maxncomp)+bad
   ovi['cf'] = dblarr(ncos,maxncomp)+bad
   ovi['tau'] = dblarr(ncos,maxncomp)+bad
   ovi['sig'] = dblarr(ncos,maxncomp)+bad
   pv = orderedhash()
   pv['ncomp'] = intarr(ncos)
   pv['weq'] = dblarr(ncos)+bad
   pv['weq_A'] = dblarr(ncos)+bad
   pv['vwtavg'] = dblarr(ncos)+bad
   pv['vwtrms'] = dblarr(ncos)+bad
   pv['v50'] = dblarr(ncos,maxncomp)+bad
   pv['cf'] = dblarr(ncos,maxncomp)+bad
   pv['tau'] = dblarr(ncos,maxncomp)+bad
   pv['sig'] = dblarr(ncos,maxncomp)+bad
   for i=0,ncos-1 do begin
      fitdir_gal=fitdir+tabdat['sgal',i]+'/'
      if file_test(fitdir_gal,/dir) then begin
         file_tmp = fitdir_gal+tabdat['sgal',i]+'NVpar_best.txt'
         if file_test(file_tmp) then begin
            readcol,file_tmp,ncomp_tmp,/silent,numline=1,skip=3,format='(I0,X)'
;            readcol,file_tmp,ncompem_tmp,/silent,numline=1,skip=1,format='(I0,X,X,X,X)'
            readcol,file_tmp,totweq_tmp,/silent,numline=1,skip=4,format='(D0,X)'
            readcol,file_tmp,vwtavg_tmp,/silent,numline=1,skip=5,format='(D0,X)'
            readcol,file_tmp,vwtrms_tmp,/silent,numline=1,skip=6,format='(D0,X)'
            nv['ncomp',i]=ncomp_tmp[0]
            if totweq_tmp[0] gt 0 then begin
               nv['weq',i]=alog10(totweq_tmp[0])
               nv['weq_A',i]=totweq_tmp[0]
            endif
            nv['vwtavg',i] = vwtavg_tmp[0]
            nv['vwtrms',i] = vwtrms_tmp[0]
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
            readcol,file_tmp,totweq_tmp,/silent,numline=1,skip=4,format='(D0,X)'
            readcol,file_tmp,vwtavg_tmp,/silent,numline=1,skip=5,format='(D0,X)'
            readcol,file_tmp,vwtrms_tmp,/silent,numline=1,skip=6,format='(D0,X)'
            ovi['ncomp',i]=ncomp_tmp[0]
            if totweq_tmp[0] gt 0 then begin
               ovi['weq',i]=alog10(totweq_tmp[0])
               ovi['weq_A',i]=totweq_tmp[0]
            endif
            ovi['vwtavg',i] = vwtavg_tmp[0]
            ovi['vwtrms',i] = vwtrms_tmp[0]
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
            readcol,file_tmp,totweq_tmp,/silent,numline=1,skip=4,format='(D0,X)'
            readcol,file_tmp,vwtavg_tmp,/silent,numline=1,skip=5,format='(D0,X)'
            readcol,file_tmp,vwtrms_tmp,/silent,numline=1,skip=6,format='(D0,X)'
            pv['ncomp',i]=ncomp_tmp[0]
            if totweq_tmp[0] gt 0 then begin
               pv['weq',i]=alog10(totweq_tmp[0])
               pv['weq_A',i]=totweq_tmp[0]
            endif
            pv['vwtavg',i] = vwtavg_tmp[0]
            pv['vwtrms',i] = vwtrms_tmp[0]
            readcol,file_tmp,cf,tau1128,lambda1128,sig,vel,/silent,$
               numline=ncomp_tmp[0],skip=10,format='(D,D,D,D,D)'
            pv['v50',i,0:ncomp_tmp[0]-1] = vel
            pv['cf',i,0:ncomp_tmp[0]-1] = cf
            pv['tau',i,0:ncomp_tmp[0]-1] = tau1128
            pv['sig',i,0:ncomp_tmp[0]-1] = sig
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
             charsize=1.75,default=4,/quiet
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
   xtit = 'log[ $\lambda$L$\down\\lambda$ (1125 $\Angstrom$) / erg s$\up-1$ ]'
   xran = [42.5d,47d]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-8000,1000]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.75,default=4,/quiet
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
      charsize=1.75,default=4,/quiet
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
      charsize=1.75,default=4,/quiet
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
      charsize=1.75,default=4,/quiet
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
      charsize=1.75,default=4,/quiet
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
             charsize=1.75,default=4,/quiet
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
             charsize=1.75,default=4,/quiet
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
             charsize=1.75,default=4,/quiet
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
   xtit = 'log[ $\lambda$L$\down\\lambda$ (1125 $\Angstrom$) / erg s$\up-1$ ]'
   xran = [42.5d,47d]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
      yran = [0.5,3.5]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)
   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.75,default=4,/quiet
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
      charsize=1.75,default=4,/quiet
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
      charsize=1.75,default=4,/quiet
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
      charsize=1.75,default=4,/quiet
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
      charsize=1.75,default=4,/quiet
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
             charsize=1.75,default=4,/quiet
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
             charsize=1.75,default=4,/quiet
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

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.2,13.2]

   ylab = 'weq'
   ytit = 'log(W$\downeq$ / $\Angstrom$)'
   yran = [-2,2]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,tabdat[xlab,inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab,iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close



   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'luv'
   xtit = 'log[ $\lambda$L$\down\\lambda$ (1125 $\Angstrom$) / erg s$\up-1$ ]'
   xran = [42.5d,47d]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',ytit=ytit,$
      position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
      color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
      position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
      yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,tabdat[xlab,inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab,iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.02]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,tabdat[xlab,inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab,iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.301,12.799]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,lagn,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,lagn,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,lagn[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,lagn[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,lmbh,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,lmbh,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,lmbh[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,lmbh[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,leddrat,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,leddrat,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,leddrat[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,leddrat[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'alphaox'
   xtit = '$\alpha$$\downox$'
   xran = [-2.2,-1.1]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,tabdat[xlab,inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab,iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lirlbol'
   xtit = 'log(L$\downIR$/L$\downbol$)'
   xran = [-0.799,0.199]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,lirlbol,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,lirlbol,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,lirlbol[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,lirlbol[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lfirlbol'
   xtit = 'log(L$\downFIR$/L$\downbol$)'
   xran = [-2,0]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,lfirlbol,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   cgoplot,lfirlbol,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,lfirlbol[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,lfirlbol[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'nhxray'
   xtit = 'log[ N(H, X-ray) / cm$\up-2$ ]'
   xran = [19d,24d]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         inz = where(x gt 0 AND x ne bad,ctnz)
         iz = where(x eq 0,ctz)
         if ctnz gt 0 then x[inz] = alog10(x[inz])+22d
         if ctz gt 0 then x[iz] = 20d
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
         if ctz gt 0 then begin
            cgoplot,x[iz]-0.07d,ynv[iz],psym=13,symsize=2,color='Red'
            cgoplot,x[iz]-0.07d,yovi[iz],psym=13,symsize=2,color='Blue'
         endif
      endif
   endfor
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shifty = (randomu(seed)+0.5d)
         shiftx = (randomu(seed)-0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shifty
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shifty
         inz = where(x gt 0 AND x ne bad,ctnz)
         iz = where(x eq 0,ctz)
         if ctnz gt 0 then x[inz] = alog10(x[inz])+22d
         if ctz gt 0 then x[iz] = 20d + shiftx
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
            if ctz gt 0 then $
               cgoplot,x[iz]-0.07d,ynv[iz],psym=13,symsize=2,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
            if ctz gt 0 then $
               cgoplot,x[iz]-0.07d,yovi[iz],psym=13,symsize=2,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'gamxray'
   xtit = '$\Gamma$ (X-ray)'
   xran = [1,3.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shifty = (randomu(seed)+0.5d)
         shiftx = (randomu(seed)-0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shifty
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shifty
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fsoftxray'
   xtit = 'log[ F(0.5-2 keV) / erg s$\up-1$ cm$\up-2$ ]'
   xran = [-14,-10]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fhardxray'
   xtit = 'log[ F(2-10 keV) / erg s$\up-1$ cm$\up-2$ ]'
   xran = [-13,-9.8]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftxray'
   xtit = 'log[ L(0.5-2 keV) / erg s$\up-1$ ]'
   xran = [42,46]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lhardxray'
   xtit = 'log[ L(2-10 keV) / erg s$\up-1$ ]'
   xran = [42.25,46.25]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'ltotxray'
   xtit = 'log[ L(0.5-10 keV) / erg s$\up-1$ ]'
   xran = [42.5,46.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftratxray'
   xtit = 'L(0.5-2 keV) / L(0.5-10 keV)'
   xran = [0.1,0.9]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lxlbol'
   xtit = 'log[ L(0.5-10 keV) / L$\downbol$ ]'
   xran = [-3.3,-0.3]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.15,0.25,0.95,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.15,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.05,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; vwtavg vs.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.2,13.2]

   ylab = 'vwtavg'
   ytit = 'Average depth-weighted velocity (km/s)'
   yran = [1000,-9000]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,tabdat[xlab,inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab,iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


   xlab = 'luv'
   xtit = 'log[ $\lambda$L$\down\\lambda$ (1125 $\Angstrom$) / erg s$\up-1$ ]'
   xran = [42.5d,47d]

   ylab = 'vwtavg'
   ytit = 'Average depth-weighted velocity (km/s)'
   yran = [1000,-9000]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,tabdat[xlab,inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab,iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'alphaox'
   xtit = '$\alpha$$\downox$'
   xran = [-2.199,-1.101]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,tabdat[xlab,inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab,iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lirlbol'
   xtit = 'log(L$\downIR$/L$\downbol$)'
   xran = [-0.799,0.199]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.75,default=4,/quiet
   cgplot,lirlbol,nv[ylab],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',ytit=ytit,$
      err_ylow=nverr,err_yhi=nverr,$
      position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,lirlbol,ovi[ylab],psym=9,symsize=1.5,$
      color='Blue',$
      err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,lirlbol[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,lirlbol[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lfirlbol'
   xtit = 'log(L$\downFIR$/L$\downbol$)'
   xran = [-2,0]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
      charsize=1.75,default=4,/quiet
   cgplot,lfirlbol,nv[ylab],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',ytit=ytit,$
      err_ylow=nverr,err_yhi=nverr,$
      position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,lfirlbol,ovi[ylab],psym=9,symsize=1.5,$
      color='Blue',$
      err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,lfirlbol[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,lfirlbol[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close



   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.05]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,tabdat[xlab,inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,tabdat[xlab,iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.301,12.799]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,lagn,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,lagn,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,lagn[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,lagn[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.501,9.499]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,lmbh,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,lmbh,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,lmbh[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,lmbh[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-1.999,0.499]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,leddrat,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   cgoplot,leddrat,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   cgoplot,leddrat[inv_nd],tmpy_nv,psym=15,symsize=1,color='Red'
   cgoplot,leddrat[iovi_nd],tmpy_ovi,psym=9,symsize=1.5,color='Blue'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'nhxray'
   xtit = 'log[ N(H, X-ray) / cm$\up-2$ ]'
   xran = [19.001d,23.999d]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickf='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         inz = where(x gt 0 AND x ne bad,ctnz)
         iz = where(x eq 0,ctz)
         if ctnz gt 0 then x[inz] = alog10(x[inz])+22d
         if ctz gt 0 then x[iz] = 20d
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
         if ctz gt 0 then begin
            cgoplot,x[iz]-0.05d,ynv[iz],psym=13,symsize=2,color='Red',$
                    err_ylow=ynverr[iz],err_yhi=ynverr[iz]
            cgoplot,x[iz]-0.05d,yovi[iz],psym=13,symsize=2,color='Blue',$
                    err_ylow=ynverr[iz],err_yhi=ynverr[iz]
         endif
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 15d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shifty = (randomu(seed)+0.5d)
         shiftx = (randomu(seed)-0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shifty
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shifty
         inz = where(x gt 0 AND x ne bad,ctnz)
         iz = where(x eq 0,ctz)
         if ctnz gt 0 then x[inz] = alog10(x[inz])+22d
         if ctz gt 0 then x[iz] = 20d + shiftx
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
            if ctz gt 0 then $
               cgoplot,x[iz]-0.07d,ynv[iz],psym=13,symsize=2,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
            if ctz gt 0 then $
               cgoplot,x[iz]-0.07d,yovi[iz],psym=13,symsize=2,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'gamxray'
   xtit = '$\Gamma$ (X-ray)'
   xran = [1,3.499]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fsoftxray'
   xtit = 'log[ F(0.5-2 keV) / erg s$\up-1$ cm$\up-2$ ]'
   xran = [-14,-10.001]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fhardxray'
   xtit = 'log[ F(2-10 keV) / erg s$\up-1$ cm$\up-2$ ]'
   xran = [-13,-9.8]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftxray'
   xtit = 'log[ L(0.5-2 keV) / erg s$\up-1$ ]'
   xran = [42,45.999]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lhardxray'
   xtit = 'log[ L(2-10 keV) / erg s$\up-1$ ]'
   xran = [42.25,46.25]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'ltotxray'
   xtit = 'log[ L(0.5-10 keV) / erg s$\up-1$ ]'
   xran = [42.5,46.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftratxray'
   xtit = 'L(0.5-2 keV) / L(0.5-10 keV)'
   xran = [0.1,0.9]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lxlbol'
   xtit = 'log[ L(0.5-10 keV) / L$\downbol$ ]'
   xran = [-3.3,-0.3]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.19,0.25,0.99,0.95],xtickform='(A1)'
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])
         ynverr = rebin([nverr[i]],tabdat['n_xray',i])
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])
         yovierr = rebin([ovierr[i]],tabdat['n_xray',i])
         cgoplot,x,ynv,psym=15,symsize=1,color='Red',err_ylow=ynverr,$
                 err_yhi=ynverr
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.3d99,1.7d99],xtit=xtit,$
          position=[0.19,0.15,0.99,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,0.09,0.19,'Undet.',chars=1.5,/norm
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         if tabdat['nvstatus',i] ne 'X' then begin
            cgoplot,x,ynv,psym=15,symsize=1,color='Red'
            cgoplot,x,ynv,/linesty,color='Red'
         endif
         if tabdat['ovistatus',i] ne 'X' then begin
            cgoplot,x,yovi,psym=9,symsize=1.5,color='Blue'
            cgoplot,x,yovi,/linesty,color='Blue'
         endif
       endif 
   endfor
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; vwtrms vs.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.2,13.2]

   ylab = 'vwtrms'
   ytit = 'Average depth-weighted velocity RMS (km/s)'
   yran = [0,3500]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgps_close



   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.05]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,tabdat[xlab],ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.3,12.8]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,lagn,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,lagn,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,lmbh,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,lmbh,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,/nomatch,$
             charsize=1.75,default=4,/quiet
   cgplot,leddrat,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,leddrat,ovi[ylab],psym=9,symsize=1.5,$
           color='Blue'
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

   amp = '&'
   dslash = '\\'
   ndat = '\nodata'
   lineofdashes = strjoin(replicate('-',62))

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
   invtotnhnz = cgsetintersection(invtot,where(tabdat['nhxray',*,0] gt 0 $
                                               AND tabdat['nhxray',*,0] ne bad),count=ctnvtotnhnz)
   iovitotnhnz = cgsetintersection(iovitot,where(tabdat['nhxray',*,0] gt 0 $
                                               AND tabdat['nhxray',*,0] ne bad),count=ctovitotnhnz)
   ibothtotnhnz = cgsetintersection(ibothtot,where(tabdat['nhxray',*,0] gt 0 $
                                               AND tabdat['nhxray',*,0] ne bad),count=ctbothtotnhnz)
   ianytotnhnz = cgsetintersection(ianytot,where(tabdat['nhxray',*,0] gt 0 $
                                               AND tabdat['nhxray',*,0] ne bad),count=ctanytotnhnz)

   invtotnhz = cgsetintersection(invtot,where(tabdat['nhxray',*,0] eq 0),count=ctnvtotnhz)
   iovitotnhz = cgsetintersection(iovitot,where(tabdat['nhxray',*,0] eq 0),count=ctovitotnhz)
   ibothtotnhz = cgsetintersection(ibothtot,where(tabdat['nhxray',*,0] eq 0),count=ctbothtotnhz)
   ianytotnhz = cgsetintersection(ianytot,where(tabdat['nhxray',*,0] eq 0),count=ctanytotnhz)

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

   invnhnz = cgsetintersection(inv,where(tabdat['nhxray',*,0] gt 0 $
                                               AND tabdat['nhxray',*,0] ne bad),count=ctnvnhnz)
   iovinhnz = cgsetintersection(iovi,where(tabdat['nhxray',*,0] gt 0 $
                                               AND tabdat['nhxray',*,0] ne bad),count=ctovinhnz)
   ibothnhnz = cgsetintersection(iboth,where(tabdat['nhxray',*,0] gt 0 $
                                               AND tabdat['nhxray',*,0] ne bad),count=ctbothnhnz)
   ianynhnz = cgsetintersection(iany,where(tabdat['nhxray',*,0] gt 0 $
                                               AND tabdat['nhxray',*,0] ne bad),count=ctanynhnz)

   invnhz = cgsetintersection(inv,where(tabdat['nhxray',*,0] eq 0),count=ctnvnhz)
   iovinhz = cgsetintersection(iovi,where(tabdat['nhxray',*,0] eq 0),count=ctovinhz)
   ibothnhz = cgsetintersection(iboth,where(tabdat['nhxray',*,0] eq 0),count=ctbothnhz)
   ianynhz = cgsetintersection(iany,where(tabdat['nhxray',*,0] eq 0),count=ctanynhz)

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

   printf,lun_tmp,'\cutinhead{$N_{\rm H}$ $>$ 0}'
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

   printf,lun_tmp,'\cutinhead{$N_{\rm H}$ $=$ 0}'
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
         printf,lun_tmp,tabdat['gal',i],amp,'N~V',amp,nv['weq_A',i],amp,$
            nv['vwtavg',i],amp,nv['vwtrms',i],amp,nv['ncomp',i],dslash,$
            format='(A12,A3,A5,A3,D6.2,A3,D9.2,A3,D8.2,A3,I3,A3)'
         labprint=1b
      endif
      if tabdat['ovistatus',i] eq 'F' OR tabdat['ovistatus',i] eq 'GF' then begin
         if labprint then labtmp = '' else labtmp = tabdat['gal',i]
         printf,lun_tmp,labtmp,amp,'O~VI',amp,ovi['weq_A',i],amp,$
            ovi['vwtavg',i],amp,ovi['vwtrms',i],amp,ovi['ncomp',i],dslash,$
            format='(A12,A3,A5,A3,D6.2,A3,D9.2,A3,D8.2,A3,I3,A3)'
         labprint=1b
      endif
      if pv['ncomp',i] ne 0 then begin
         if labprint then labtmp = '' else labtmp = tabdat['gal',i]
         printf,lun_tmp,labtmp,amp,'P~V',amp,pv['weq_A',i],amp,$
            pv['vwtavg',i],amp,pv['vwtrms',i],amp,pv['ncomp',i],dslash,$
            format='(A12,A3,A5,A3,D6.2,A3,D9.2,A3,D8.2,A3,I3,A3)'
      endif
      
   endfor
   
   free_lun,lun_tmp

   ; Fit results

   openw,lun_tmp,tabdir+'tab_table1quants.tex',/get_lun

   printf,lun_tmp,'#Col 2: log(Lbol/Lsun) [corrected to cosmology for paper]
   printf,lun_tmp,'#Col 3: log(LIR/Lsun) [corrected to cosmology for paper]
   printf,lun_tmp,'#Col 4: log(lam L_lam) at 1125 A; 1.5% bandpass; corrected for Gal. extinction
   printf,lun_tmp,'#Col 5: log(MBH/Msun)
   printf,lun_tmp,'#Col 6: log(Edd.rat.)
   for i=0,ncos-1 do begin
      printf,lun_tmp,tabdat['gal',i],tabdat['lbol',i],tabdat['lir',i],$
         tabdat['luv',i],lmbh[i],leddrat[i],format='(A-12,D8.2,D8.2,D8.2,D8.2,D8.2)'
   endfor

   free_lun,lun_tmp

end
