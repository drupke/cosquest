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
;    Copyright (C) 2016 David S. N. Rupke
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

   qsotab='/Users/drupke/Dropbox/qsos/qsos.csv'
   cosdir='/Users/drupke/Box Sync/cosquest/'
   fitdir=cosdir+'fits/'
   plotdir=cosdir+'plots/correlations/'

;  Read table
   rows=[3,81]
   gal = read_csvcol(qsotab,'A',rows=rows,sep=',',type='string')
   sgal = read_csvcol(qsotab,'C',rows=rows,sep=',',type='string')
   z = read_csvcol(qsotab,'D',rows=rows,sep=',',junk=bad)
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

;  Parse table data
   icos = where(cossamp eq 'V18',ncos)
   nlines = n_elements(gal)
   tabdat = hash()
   tabdat['gal'] = gal[icos]
   tabdat['sgal'] = sgal[icos]
   tabdat['z'] = z[icos]
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
            tabdat['fsoftxray',j,k] = fsoftxray[i]
            tabdat['fsoftxray_errlo',j,k] = fsoftxray_errlo[i]
            tabdat['fsoftxray_errhi',j,k] = fsoftxray_errhi[i]
            tabdat['fhardxray',j,k] = fhardxray[i]
            tabdat['fhardxray_errlo',j,k] = fhardxray_errlo[i]
            tabdat['fhardxray_errhi',j,k] = fhardxray_errhi[i]
            tabdat['lsoftxray',j,k] = lsoftxray[i]
            tabdat['lhardxray',j,k] = lhardxray[i]
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
               tabdat['fsoftxray',j,k] = fsoftxray[i]
               tabdat['fsoftxray_errlo',j,k] = fsoftxray_errlo[i]
               tabdat['fsoftxray_errhi',j,k] = fsoftxray_errhi[i]
               tabdat['fhardxray',j,k] = fhardxray[i]
               tabdat['fhardxray_errlo',j,k] = fhardxray_errlo[i]
               tabdat['fhardxray_errhi',j,k] = fhardxray_errhi[i]
               tabdat['lsoftxray',j,k] = lsoftxray[i]
               tabdat['lhardxray',j,k] = lhardxray[i]
               k++
            endif else begin
               k = 0 ; re-zero x-ray component index
            endelse
         endif
      endelse
   endfor

; Compute physical quantities
; 
;  Set AGN fraction to 1 if we don't have a measurement
   ibdagnfrac = where(tabdat['agnfrac'] eq bad,ctbdagnfrac)
   if ctbdagnfrac gt 0 then begin
      tabdat['agnfrac',ibdagnfrac] = 1d
      tabdat['agnfraclb',ibdagnfrac] = 1d
      tabdat['agnfracub',ibdagnfrac] = 1d
   endif

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


;  Get components
   maxncomp = 15
   nv = hash()
   nv['ncomp'] = intarr(ncos)
   nv['weq'] = dblarr(ncos)+bad
   nv['vwtavg'] = dblarr(ncos)+bad
   nv['vwtrms'] = dblarr(ncos)+bad
   nv['v50'] = dblarr(ncos,maxncomp)+bad
   nv['cf'] = dblarr(ncos,maxncomp)+bad
   nv['tau'] = dblarr(ncos,maxncomp)+bad
   nv['sig'] = dblarr(ncos,maxncomp)+bad
   ovi = hash()
   ovi['ncomp'] = intarr(ncos)
   ovi['weq'] = dblarr(ncos)+bad
   ovi['vwtavg'] = dblarr(ncos)+bad
   ovi['vwtrms'] = dblarr(ncos)+bad
   ovi['v50'] = dblarr(ncos,maxncomp)+bad
   ovi['cf'] = dblarr(ncos,maxncomp)+bad
   ovi['tau'] = dblarr(ncos,maxncomp)+bad
   ovi['sig'] = dblarr(ncos,maxncomp)+bad
   for i=0,ncos-1 do begin
      fitdir_gal=fitdir+tabdat['sgal',i]+'/'
      if file_test(fitdir_gal,/dir) then begin
         file_tmp = fitdir_gal+tabdat['sgal',i]+'NVpar_best.txt'
         if file_test(file_tmp) then begin
            readcol,file_tmp,ncomp_tmp,/silent,numline=1,format='(I0,X,X,X,X)'
            readcol,file_tmp,ncompem_tmp,/silent,numline=1,skip=1,format='(I0,X,X,X,X)'
            readcol,file_tmp,totweq_tmp,/silent,numline=1,skip=2,format='(D0,X)'
            readcol,file_tmp,vwtavg_tmp,/silent,numline=1,skip=3,format='(D0,X)'
            readcol,file_tmp,vwtrms_tmp,/silent,numline=1,skip=4,format='(D0,X)'
            nv['ncomp',i]=ncomp_tmp[0]
            if totweq_tmp[0] gt 0 then nv['weq',i]=alog10(totweq_tmp[0])
            nv['vwtavg',i] = vwtavg_tmp[0]
            nv['vwtrms',i] = vwtrms_tmp[0]
            readcol,file_tmp,cf,tau1243,lambda1243,sig,vel,/silent,$
                    numline=ncomp_tmp[0],skip=6,format='(D,D,D,D,D)'
            nv['v50',i,0:ncomp_tmp[0]-1] = vel
            nv['cf',i,0:ncomp_tmp[0]-1] = cf
            nv['tau',i,0:ncomp_tmp[0]-1] = tau1243
            nv['sig',i,0:ncomp_tmp[0]-1] = sig
         endif
         file_tmp = fitdir_gal+tabdat['sgal',i]+'OVIpar_best.txt'
         if file_test(file_tmp) then begin
            readcol,file_tmp,ncomp_tmp,/silent,numline=1,format='(I0,X,X,X,X)'
            readcol,file_tmp,ncompem_tmp,/silent,numline=1,skip=1,format='(I0,X,X,X,X)'
            readcol,file_tmp,totweq_tmp,/silent,numline=1,skip=2,format='(D0,X)'
            readcol,file_tmp,vwtavg_tmp,/silent,numline=1,skip=3,format='(D0,X)'
            readcol,file_tmp,vwtrms_tmp,/silent,numline=1,skip=4,format='(D0,X)'
            ovi['ncomp',i]=ncomp_tmp[0]
            if totweq_tmp[0] gt 0 then ovi['weq',i]=alog10(totweq_tmp[0])
            ovi['vwtavg',i] = vwtavg_tmp[0]
            ovi['vwtrms',i] = vwtrms_tmp[0]
            readcol,file_tmp,cf,tau1038,lambda1038,sig,vel,/silent,$
                    numline=ncomp_tmp[0],skip=6,format='(D,D,D,D,D)'
            ovi['v50',i,0:ncomp_tmp[0]-1] = vel
            ovi['cf',i,0:ncomp_tmp[0]-1] = cf
            ovi['tau',i,0:ncomp_tmp[0]-1] = tau1038
            ovi['sig',i,0:ncomp_tmp[0]-1] = sig
         endif
      endif
   endfor
   igd_nv_comp = where(nv['v50'] ne bad,ctnvcomp)
   igd_ovi_comp = where(ovi['v50'] ne bad,ctovicomp)



;  Plots


tvlct,[[27],[158],[119]],100
tvlct,[[217],[95],[2]],101
tvlct,[[117],[112],[179]],102
tvlct,[[231],[41],[138]],103
tvlct,[[102],[166],[30]],104
tvlct,[[230],[171],[2]],105
tvlct,[[166],[118],[29]],106
tvlct,[[102],[102],[102]],107
tvlct,[[0],[0],[0]],108

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'z'
   xtit = 'Redshift'
   xran = [0d,0.4d]

   ylab = 'lbol'
   ytit = 'log(L$\downbol$/L$\sun$)'
   yran = [11.3,12.8]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,tabdat[xlab],tabdat[ylab],xran=xran,yran=yran,$
          psym=16,symsize=1,color='Black',xtit=xtit,ytit=ytit

   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-8000,1000]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)
   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=16,symsize=1,$
           color='Blue'
;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close



   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.05]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-8000,1000]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)
   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=16,symsize=1,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-8000,1000]

   xrebin = rebin(lagn,ncos,maxncomp)
   xrebin = rebin(lagn,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=16,symsize=1,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.5]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-8000,1000]

   xrebin = rebin(lmbh,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=16,symsize=1,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.5]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-8000,1000]

   xrebin = rebin(leddrat,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=16,symsize=1,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'cf'
   xtit = 'Covering Factor'
   xran = [-0.1,1.1]

   ylab = 'v50'
   ytit = 'v$\down50$ (km/s)'
   yran = [-9000,1000]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,nv[xlab,igd_nv_comp],nv[ylab,igd_nv_comp],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,ovi[xlab,igd_ovi_comp],ovi[ylab,igd_ovi_comp],psym=16,symsize=1,$
           color='Blue'
;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
   yran = [0.5,3.5]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)
   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=16,symsize=1,$
           color='Blue'
;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close



   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.05]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
   yran = [0.5,3.5]

   xrebin = rebin(tabdat[xlab],ncos,maxncomp)
   xrebin = rebin(tabdat[xlab],ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=16,symsize=1,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
   yran = [0.5,3.5]

   xrebin = rebin(lagn,ncos,maxncomp)
   xrebin = rebin(lagn,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=16,symsize=1,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.5]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
   yran = [0.5,3.5]

   xrebin = rebin(lmbh,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=16,symsize=1,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.5]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
   yran = [0.5,3.5]

   xrebin = rebin(leddrat,ncos,maxncomp)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1,default=4,/quiet
   cgplot,xrebin[igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
      psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   ;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
   ;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,xrebin[igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=16,symsize=1,$
      color='Blue'
   ;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
   ;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close
   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'cf'
   xtit = 'Covering Factor'
   xran = [-0.1,1.1]

   ylab = 'sig'
   ytit = 'log[ $\sigma$ / km/s ]'
   yran = [0.5,3.5]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,nv[xlab,igd_nv_comp],alog10(nv[ylab,igd_nv_comp]),xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
;          err_ylo=nv['sig',igd_nv_comp]*fwhm2sig/2d,$
;          err_yhi=nv['sig',igd_nv_comp]*fwhm2sig/2d
   cgoplot,ovi[xlab,igd_ovi_comp],alog10(ovi[ylab,igd_ovi_comp]),psym=16,symsize=1,$
           color='Blue'
;           err_ylo=ovi['sig',igd_ovi_comp]*fwhm2sig/2d,$
;           err_yhi=ovi['sig',igd_ovi_comp]*fwhm2sig/2d
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'weq'
   ytit = 'log(W$\downeq$ / A)'
   yran = [-2,2]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,tabdat[xlab],ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close



   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.05]

   ylab = 'weq'
   ytit = 'log(W$\downeq$ / A)'
   yran = [-2,2]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,tabdat[xlab],ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'weq'
   ytit = 'log(W$\downeq$ / A)'
   yran = [-2,2]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,lagn,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,lagn,ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.5]

   ylab = 'weq'
   ytit = 'log(W$\downeq$ / A)'
   yran = [-2,2]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,lmbh,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,lmbh,ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.5]

   ylab = 'weq'
   ytit = 'log(W$\downeq$ / A)'
   yran = [-2,2]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,leddrat,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,leddrat,ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'vwtavg'
   ytit = 'Average depth-weighted velocity (km/s)'
   yran = [1000,-9000]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr
   cgoplot,tabdat[xlab],ovi[ylab],psym=16,symsize=1,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.05]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr
   cgoplot,tabdat[xlab],ovi[ylab],psym=16,symsize=1,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.3,12.8]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,lagn,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr
   cgoplot,lagn,ovi[ylab],psym=16,symsize=1,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.5]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,lmbh,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr
   cgoplot,lmbh,ovi[ylab],psym=16,symsize=1,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.5]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,leddrat,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit,$
          err_ylow=nverr,err_yhi=nverr
   cgoplot,leddrat,ovi[ylab],psym=16,symsize=1,$
           color='Blue',$
           err_ylow=ovierr,err_yhi=ovierr
   cgoplot,xran,[0,0],/linesty
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'nhxray'
   xtit = 'log[ N(H, X-ray) / cm$\up-2$ ]'
   xran = [19.5d,24d]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,xtit=xtit,ytit=ytit
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
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
         if ctz gt 0 then begin
            cgoplot,x[iz]-0.05d,ynv[iz],psym=20,symsize=2,color='Red',$
                    err_ylow=ynverr[iz],err_yhi=ynverr[iz]
            cgoplot,x[iz]-0.05d,yovi[iz],psym=20,symsize=2,color='Blue',$
                    err_ylow=ynverr[iz],err_yhi=ynverr[iz]
         endif
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'gamxray'
   xtit = '$\Gamma$ (X-ray)'
   xran = [1,3.5]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.12,0.25,0.95,0.95],xtickform='(A1)'
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
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.1d99,1.9d99],xtit=xtit,$
          position=[0.12,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1
   cgtext,1.06,1.4d99,'No UV lines'
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
       endif 
   endfor
   cgps_close



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fsoftxray'
   xtit = 'log[ F(0.5-2 keV) / erg s$\up-1$ cm$\up-2$ ]'
   xran = [-14,-8]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.12,0.25,0.95,0.95],xtickform='(A1)'
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
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.1d99,1.9d99],xtit=xtit,$
          position=[0.12,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,-13.9,1.4d99,'No UV lines'
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])-12d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'fhardxray'
   xtit = 'log[ F(2-10 keV) / erg s$\up-1$ cm$\up-2$ ]'
   xran = [-13,-8]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.12,0.25,0.95,0.95],xtickform='(A1)'
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
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.1d99,1.9d99],xtit=xtit,$
          position=[0.12,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,-13.9,1.4d99,'No UV lines'
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]-12d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lsoftxray'
   xtit = 'log[ L(0.5-2 keV) / erg s$\up-1$ ]'
   xran = [42,48]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.12,0.25,0.95,0.95],xtickform='(A1)'
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
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.1d99,1.9d99],xtit=xtit,$
          position=[0.12,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,42.2,1.4d99,'No UV lines'
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = alog10(tabdat[xlab,i,0:tabdat['n_xray',i]-1])+44d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lhardxray'
   xtit = 'log[ L(2-10 keV) / erg s$\up-1$ ]'
   xran = [43,46]

   nverr = nv['vwtrms']
   ibdnverr = where(nverr eq bad,ctbd)
   if ctbd gt 0 then nverr[ibdnverr] = 0d
   ovierr = ovi['vwtrms']
   ibdovierr = where(ovierr eq bad,ctbd)
   if ctbd gt 0 then ovierr[ibdovierr] = 0d

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,[0],/nodat,xran=xran,yran=yran,ytit=ytit,$
          position=[0.12,0.25,0.95,0.95],xtickform='(A1)'
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
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue',err_ylow=yovierr,$
                 err_yhi=yovierr
         cgoplot,x,yovi,/linesty,color='Blue'
      endif
   endfor
   cgoplot,xran,[0,0],/linesty
   cgplot,[0],/nodat,xran=xran,yran=[0.1d99,1.9d99],xtit=xtit,$
          position=[0.12,0.15,0.95,0.25],/noerase,/xsty,/ysty,$
          yticks=1,ytickf='(A1)'
   cgtext,42.2,1.4d99,'No UV lines'
   seed = 5d
   for i=0,ncos-1 do begin
      if tabdat['n_xray',i] ge 1 then begin
         x = tabdat[xlab,i,0:tabdat['n_xray',i]-1]+44d
         shift = (randomu(seed)+0.5d)
         ynv = rebin([nv[ylab,i]],tabdat['n_xray',i])*shift
         yovi = rebin([ovi[ylab,i]],tabdat['n_xray',i])*shift
         cgoplot,x,ynv,psym=15,symsize=1,color='Red'
         cgoplot,x,ynv,/linesty,color='Red'
         cgoplot,x,yovi,psym=16,symsize=1,color='Blue'
         cgoplot,x,yovi,/linesty,color='Blue'
       endif 
   endfor
   cgps_close


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lbol'
   xtit = 'log(L$\downbol$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'vwtrms'
   ytit = 'Average depth-weighted velocity RMS (km/s)'
   yran = [0,3500]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,tabdat[xlab],ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close



   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'agnfrac'
   xtit = 'AGN fraction'
   xran = [0.7,1.05]

   ylab = 'vwtrms'
   ytit = 'Average depth-weighted velocity RMS (km/s)'
   yran = [0,3500]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,tabdat[xlab],nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,tabdat[xlab],ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lagn'
   xtit = 'log(L$\downAGN$/L$\sun$)'
   xran = [11.3,12.8]

   ylab = 'vwtrms'
   ytit = 'Average depth-weighted velocity RMS (km/s)'
   yran = [0,3500]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,lagn,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,lagn,ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'lmbh'
   xtit = 'log(M$\downBH$/M$\sun$)'
   xran = [6.5,9.5]

   ylab = 'vwtrms'
   ytit = 'Average depth-weighted velocity RMS (km/s)'
   yran = [0,3500]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,lmbh,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,lmbh,ovi[ylab],psym=16,symsize=1,$
           color='Blue'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   xlab = 'leddrat'
   xtit = 'Eddington Ratio'
   xran = [-2,0.5]

   ylab = 'vwtrms'
   ytit = 'Average depth-weighted velocity RMS (km/s)'
   yran = [0,3500]

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
             charsize=1,default=4,/quiet
   cgplot,leddrat,nv[ylab],xran=xran,yran=yran,$
          psym=15,symsize=1,color='Red',xtit=xtit,ytit=ytit
   cgoplot,leddrat,ovi[ylab],psym=16,symsize=1,$
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
   ig_nv = where(nvstatus_sort eq 'G',ct_g_nv)
   ig_ovi = where(ovistatus_sort eq 'G',ct_g_ovi)
   
   igd_nv_comp_sort = where(nv_v50_sort ne bad)
   igd_ovi_comp_sort = where(ovi_v50_sort ne bad)
   sigthresh = 25d
   igd_nv_comp_sort_nar = where(nv_sig_sort ne bad AND nv_sig_sort le sigthresh)
   igd_ovi_comp_sort_nar = where(ovi_sig_sort ne bad AND ovi_sig_sort le sigthresh)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=10,$
             charsize=1,default=4,/quiet

   cgplot,nv_v50_sort[igd_nv_comp_sort],yrebin[igd_nv_comp_sort],$
          xran=xran,yran=yran,xtit=xtit,ytit=ytit,$
          psym=15,symsize=1,color='Red',$
          ytickname=gal_sort,ytickv=indgen(ncos)+1,yticks=ncos-1,$
          position=[0.2,0.15,0.8,0.99],/noerase
   cgoplot,nv_v50_sort[igd_nv_comp_sort_nar],yrebin[igd_nv_comp_sort_nar],$
          psym=6,symsize=2,color='Red'
   cgoplot,ovi_v50_sort[igd_ovi_comp_sort],yrebin[igd_ovi_comp_sort],$
           psym=16,symsize=1,color='Blue'
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
          ytickf='(A1)',xticks=3,xtickv=[1,2],xtickname=['[NV]','[OVI]']
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

;   al_legend,['[N V] 1238,1243','[O VI] 1032,1038','not in spectral range',$
;      'affected by geocoronal line or chip gap','undetected'],$
;      /norm,spacing=1.5,$
;      color=['Red','Blue','Black','Black','Black'],$
;      psym=[6,16,7,36,11],symsize=[1.5,1,1,1.5,1],position=[0.555,0.12]
   al_legend,['[N V] 1238,1243','[O VI] 1032,1038',$
              '[N V], $\sigma$ < 25 km/s',$
              '[O VI], $\sigma$ < 25 km/s'],$
              /norm,spacing=1.5,$
              color=['Red','Blue','Red','Blue'],$
              psym=[15,16,6,9],symsize=[1,1,2,2],position=[0.2,0.1],margin=.6
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
   ig_nv = where(nvstatus_sort eq 'G',ct_g_nv)
   ig_ovi = where(ovistatus_sort eq 'G',ct_g_ovi)
   
   igd_nv_comp_sort = where(nv_v50_sort ne bad)
   igd_ovi_comp_sort = where(ovi_v50_sort ne bad)
   sigthresh = 25d
   igd_nv_comp_sort_nar = where(nv_sig_sort ne bad AND nv_sig_sort le sigthresh)
   igd_ovi_comp_sort_nar = where(ovi_sig_sort ne bad AND ovi_sig_sort le sigthresh)

   cgps_open,plotdir+ylab+'_vs_'+xlab+'.eps',/encap,/inches,xsiz=7.5,ysize=10,$
             charsize=1,default=4,/quiet

   cgplot,nv_v50_sort[igd_nv_comp_sort],yrebin[igd_nv_comp_sort],$
          xran=xran,yran=yran,xtit=xtit,ytit=ytit,$
          psym=15,symsize=1,color='Red',$
          ytickname=gal_sort,ytickv=indgen(ncos)+1,yticks=ncos-1,$
          position=[0.2,0.15,0.8,0.99],/noerase
   cgoplot,nv_v50_sort[igd_nv_comp_sort_nar],yrebin[igd_nv_comp_sort_nar],$
          psym=6,symsize=2,color='Red'
   cgoplot,ovi_v50_sort[igd_ovi_comp_sort],yrebin[igd_ovi_comp_sort],$
           psym=16,symsize=1,color='Blue'
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
          ytickf='(A1)',xticks=3,xtickv=[1,2],xtickname=['[NV]','[OVI]']
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

;   al_legend,['[N V] 1238,1243','[O VI] 1032,1038','not in spectral range',$
;      'affected by geocoronal line or chip gap','undetected'],$
;      /norm,spacing=1.5,$
;      color=['Red','Blue','Black','Black','Black'],$
;      psym=[6,16,7,36,11],symsize=[1.5,1,1,1.5,1],position=[0.555,0.12]
   al_legend,['[N V] 1238,1243','[O VI] 1032,1038',$
              '[N V], $\sigma$ < 25 km/s',$
              '[O VI], $\sigma$ < 25 km/s'],$
              /norm,spacing=1.5,$
              color=['Red','Blue','Red','Blue'],$
              psym=[15,16,6,9],symsize=[1,1,2,2],position=[0.2,0.1],margin=.6
   al_legend,['not in spectral range',$
              'affected by geocoronal line or chip gap','undetected'],$
              /norm,spacing=1.5,$
              color=['Black','Black','Black'],$
              psym=[7,36,11],symsize=[1,1.5,1],position=[0.5,0.1],margin=.6

   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'hist_cf.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
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
  
   cgplot,[0],xran=[-0.1,1.1],yran=[0,0.5],/ysty,xticks=4,xtickv=[0.1,0.3,0.5,0.7,0.9],$
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



   al_legend,['NV','OIV','NV+OVI'],$
              color=[101,102,103],psym=[15,15,15],$
              position=[-0.1,0.5],spacing=2,/data,symsize=[2,2,2],back='White'
   cgps_close

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'hist_sig.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
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

   al_legend,['NV','OIV','NV+OVI'],$
              color=[101,102,103],psym=[15,15,15],$
              position=[0,0.8],spacing=2,/data,symsize=[2,2,2],back='White'
   cgps_close

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'hist_tau.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1.5,default=4,/quiet

   bin1 = cghistogram(alog10(nv['tau',igd_nv_comp]),binsize=0.6d,min=-1.6d,missing=bad,/freq)
   bin2 = cghistogram(alog10(ovi['tau',igd_ovi_comp]),binsize=0.6d,min=-1.6d,missing=bad,/freq)
   bin3 = (bin1+bin2)/2d
   bin1n = cghistogram(alog10(nv['tau',igd_nv_comp]),binsize=0.6d,min=-1.6d,missing=bad)
   bin2n = cghistogram(alog10(ovi['tau',igd_ovi_comp]),binsize=0.6d,min=-1.6d,missing=bad)
   bin3n = bin1n+bin2n

   cgplot,[0],xran=[0,5],yran=[0,0.8],/ysty,xticks=3,xtickv=[1,2,3,4],$
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

   al_legend,['NV','OIV','NV+OVI'],$
      color=[101,102,103],psym=[15,15,15],$
      position=[0,0.8],spacing=2,/data,symsize=[2,2,2],back='White'
   cgps_close


   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   cgps_open,plotdir+'hist_vel.eps',/encap,/inches,xsiz=7.5,ysize=7.5,$
      charsize=1.5,default=4,/quiet

   bin1 = cghistogram(nv['v50',igd_nv_comp],binsize=2000d,min=-8000d,missing=bad,/freq)
   bin2 = cghistogram(ovi['v50',igd_ovi_comp],binsize=2000d,min=-8000d,missing=bad,/freq)
   bin3 = (bin1+bin2)/2d
   bin1n = cghistogram(nv['v50',igd_nv_comp],binsize=2000d,min=-8000d,missing=bad)
   bin2n = cghistogram(ovi['v50',igd_ovi_comp],binsize=2000d,min=-8000d,missing=bad)
   bin3n = bin1n+bin2n

   cgplot,[0],xran=[0,6],yran=[0,0.8],/ysty,xticks=4,xtickv=[1,2,3,4,5],$
      xtickname=['-8000 to -6000','-6000 to -4000','-4000 to -2000',$
                 '-2000 to 0','> 0',' '],$
      /nodata,xminor=1,title='Velocity (km/s) Distributions',xcharsize=0.5

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

   al_legend,['NV','OIV','NV+OVI'],$
      color=[101,102,103],psym=[15,15,15],$
      position=[0,0.8],spacing=2,/data,symsize=[2,2,2],back='White'
   cgps_close

end
