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
pro cos_corr

   bad = 1d99
   fwhm2b = 2d*sqrt(alog(2d))
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
   loglbol = read_csvcol(qsotab,'G',rows=rows,sep=',',junk=bad)
   agnfrac = read_csvcol(qsotab,'H',rows=rows,sep=',',junk=bad)
   mbh = read_csvcol(qsotab,'K',rows=rows,sep=',',junk=bad)
   cossamp = read_csvcol(qsotab,'Y',rows=rows,sep=',',type='string')
   nhxray = read_csvcol(qsotab,'AM',rows=rows,sep=',',junk=bad)

;  Parse table data
   icos = where(cossamp eq 'V16',ncos)
   nlines = n_elements(gal)
   tabdat = hash()
   tabdat.gal = gal[icos]
   tabdat.sgal = sgal[icos]
   tabdat.z = z[icos]
   tabdat.loglbol = loglbol[icos]
;   tabdat.nxray = dblarr(ngal)
;   tabdat.nhxray = dblarr(ngal,6)+bad
;   j = -1 ; galaxy index
;   for i=0,gal-1 do begin
;      if gal[i] ne '' then begin
;         j++ ; move to next galaxy
;         tabdat.nxray = k+1 ; record number of x-ray components
;         k=0 ; zero x-ray component index for next galaxy
;         tabdat.nhxray[j,k] = nhxray[i] ; record first x-ray measurement
;      endif else begin
;         k++
;         tabdat.nhxray[
;      endelse
;   endfor

;  Get components
   nv = hash()
   nv.ncomp = intarr(ncos)
   nv.v50 = dblarr(ncos,15)+bad
   ovi = hash()
   ovi.ncomp = intarr(ncos)
   ovi.v50 = dblarr(ncos,15)+bad
   for i=0,ncos-1 do begin
      fitdir_gal=fitdir+tabdat.sgal[i]+'/'
      if file_exist(fitdir,/dir) then begin
         file_tmp = fitdir+tabdat.sgal[i]+'NVpar_best.txt'
         if file_exist(file_tmp) then begin
            readcol,file_tmp,ncomp_tmp,/silent,numline=1,format='(I0,X,X,X,X)'
            nv.ncomp[i]=ncomp_tmp[0]
            readcol,file_tmp,cf,tau1243,lambda1243,
         endif
      endif
   endfor

   x = 'nv_v50'
   y = 'lbol'
   cgps_open,plotdir+x+'_vs_'+y+'.eps',/encap,/inches,xsiz=7.5,ysize=7.5
   cgps_close

end