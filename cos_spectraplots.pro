; docformat = 'rst'
;
;+
;
; Produces plots showing the entire data spectra, as well as zoomed-in
; portions.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Postscript file with plots.
;
; :Params:
;    table: in, required, type=str
;      File that holds a list of galaxies and redshifts.
;    infile: in, required, type=str
;      Input spectrum as a text file.
;    outfile: in, required, type=str
;      Filename of output plot.
;    galaxyshortname: in, required, type=str
;      Galaxy being plotted.
;
; :Keywords:
;    ignorewave: in, optional, type=dblarr(N,2)
;      List of N regions to ignore in determining vertical plot range. Each 
;      region has a lower and upper wavelength limit.
;
; :Author:
;    Anthony Dinh To::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104  
;      andto94@gmail.com
;
; :History:
;    ChangeHistory::
;      2016jul06, ADT, created
;      2016jul22, DSNR, cosmetic and input changes
;      2016sep07, DSNR, changed line label logic; moved line list to 
;                       IFSF_LINELIST; added logic to prevent label collisions
;      2020dec18, DSNR, fix for new IFSF_LINELIST logic; added S/N calculation
;                       and output
;
; :Copyright:
;    Copyright (C) 2016--2020 Anthony To, David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
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
PRO cos_spectraplots, table, infile, outfile, galaxyshortname,$
                      ignorewave=ignorewave,twave=twave,tflux=tflux,$
                      medsn=medsn,flam=flam

;  List of emission/absorption lines and their corresponding wavelengths.
   linelab = 1b
   lines = ifsf_linelist(!NULL,linelab=linelab,/vacuum,/all,/quiet)
;  Get lists from hashes. Make sure that values and labels line up.
   LineWavelength = (lines.values()).toarray()
   LineLabel = strarr(n_elements(lines))
   for i=0,n_elements(lines)-1 do $
      LineLabel[i] = linelab[(lines.keys())[i]]

; Read wavelength and flux
  readcol, infile, wave, flux, err, /silent
  
  if keyword_set(medsn) then begin
  ; Get S/N near 1300
    llow = value_locate(wave,1290d)
    lhi = value_locate(wave,1310d)
    medsn = median(flux[llow:lhi]/err[llow:lhi])
    if llow ne -1 and lhi ne -1 then $
      print,'Median S/N in range 1290-1310 A: ',medsn,format='(A0,D0.2)'
  endif

; Read galaxy full names and redshifts
  trows=[3,85]
  name = read_csvcol(table,'A',rows=trows,sep=',',type='string')
  galaxyshortnamelist = read_csvcol(table,'C',rows=trows,sep=',',type='string')
  z = read_csvcol(table,'D',rows=trows,sep=',',junk=bad)

  selectionparameter=WHERE(galaxyshortnamelist eq galaxyshortname)
  galaxyfullname=name[selectionparameter[0]]
  zsys=z[selectionparameter[0]]

  if keyword_set(flam) then begin
     ; Get flux at rest-frame 1125 A
     lam1125obs = 1125*(1d + zsys)
     bandpass = 0.015d * lam1125obs
     llow = value_locate(wave,lam1125obs - bandpass/2d)
     lhi = value_locate(wave,lam1125obs + bandpass/2d)
     flam = mean(flux[llow:lhi])
  endif

  
; Shifts the absorption/emission waveelngths by the galaxy's redshift
  ShiftedLines= LineWavelength*(1+zsys)
  
; Avoid line label collisions
  closethresh = 0.3d
  nlines = n_elements(LineWavelength)
  linesall = [LineWavelength,ShiftedLines]
  linelaball = [LineLabel,LineLabel]
  linecolall = [strarr(nlines)+'Blue',strarr(nlines)+'Red']
  lineyfracall = [dblarr(nlines)+0.05d,dblarr(nlines)+0.05d]
  isort_linesall = sort(linesall)
  sort_linesall = linesall[isort_linesall]
  sort_linelaball = linelaball[isort_linesall]
  sort_linecolall = linecolall[isort_linesall]
  sort_lineyfracall = lineyfracall[isort_linesall]
  dlines = sort_linesall[1:nlines*2-1] - sort_linesall[0:nlines*2-2]
  iclose = where(dlines lt closethresh,ctclose)
  if ctclose gt 0 then begin
     for i=0,ctclose-1 do begin
         sort_linelaball[iclose[i]] = $
            sort_linelaball[iclose[i]]+', '+sort_linelaball[iclose[i]+1]
         sort_linelaball[iclose[i]+1] = ''
         sort_lineyfracall[iclose[i]] = 0.05d
     endfor
  endif
  
; Rescaling flux values
  relativeflux=flux/1E-14

; Get rid of geocoronal Ly-alpha when scaling plots
  igcly_not = where(wave le 1212d OR wave ge 1220d)

; Get rid of other regions
  if keyword_set(ignorewave) then begin
     signore = size(ignorewave)
     if signore[0] eq 1 then nignore=1 else nignore=signore[2]
     for i=0,nignore-1 do begin
        igd_tmp = where(wave le ignorewave[0,i] OR wave ge ignorewave[1,i])
        if i eq 0 then iignore_not = igd_tmp else $
        iignore_not = cgsetintersection(iignore_not,igd_tmp)
     endfor
     igdwave = cgsetintersection(iignore_not,igcly_not)
  endif else igdwave = igcly_not

; Arrays with only good regions
  relativeflux_gd = relativeflux[igdwave]
  wave_gd = wave[igdwave]
  
; Set the y-range of the big spectra
  yran=[0,Max(relativeflux_gd)]

; Acquires the range of the wavelength values
  baserange=Max(wave)-Min(wave)

  xsize=7.5
  ysize=8.5 
  xoffset=(8.5-xsize)/2.0d
  yoffset=(11.0-ysize)/2.0d

; Various plotting defaults  
  plotchars=0.75
  labchars=0.5
  labchart=1
  x1=0.07
  x2=0.99
  
;Opens a Postscript value to draw plots on
  cgPS_OPEN,outfile,$
            /Encapsulated,scale_factor=1,Charsize=.5,/NOMATCH,$
            xsize=xsize, ysize=ysize, xoffset=xoffset, yoffset=yoffset, /Inches
  !Y.THICK=2
  !X.THICK=2
  !Y.TICKFORMAT='(F0.1)'
  !Y.MINOR = 2
  !X.TICKLEN=0.05
  !Y.TICKLEN=0.01
    
  cgText,.5,.98,galaxyfullname+', '+'z = '+String(zsys,format='(D0.3)'), $
         alignment=.5, Charsize = 1,/norm
  
;Full range plot
  cgplot, wave, relativeflux, xstyle=1, ystyle=1, yran=yran,$
          xtit='Observed Wavelength ($\Angstrom$)',$
          ytit='F!I$\lambda$!N/10!E-14!N (ergs s!E-1 !Ncm!E-2 !N$\Angstrom$!E-1!N)', $
          Position = [x1,.86,x2,.97], CHARSIZE=plotchars,thick=0.5,/NoErase,$
          xticklen=0.05,yticklen=0.01

  if keyword_set(twave) AND keyword_set(tflux) then $
     cgoplot,twave,tflux,color='Red'
    
;Plot legend
  AL_LEGEND,['Intrinsic','Galactic'], $
            Color=['Red','Blue'], charsize=1, charthick=2, $
            Linestyle=[0,0], $
            Position = [Min(wave)+.04*(Max(wave)-Min(wave)),$
                        .9*(Max(relativeflux)-Min(relativeflux))], $
            bthick=.6,clear=1,linsize=0.25


; Plots 6 zoomed-in regions, one after another, $
; along with absorption/emission lines and labels

   y1 = [0.695,0.56,0.425,0.29,0.155,0.02]
   y2 = [0.805,0.67,0.535,0.40,0.265,0.13]
   for i=0,5 do begin
      xran = [Min(wave)+(double(i)/6)*baserange,$
              Min(wave)+(double(i+1)/6)*baserange]
      yran = [0,3*Sqrt(Mean((relativeflux_gd[value_locate(wave_gd,xran[0]):$
                                             value_locate(wave_gd,xran[1])])^2d))]
      cgplot, wave, relativeflux, xstyle=1, ystyle=1, $
              xran=xran, axiscolor='Black',color='Black',yran=yran, $
              Position = [x1,y1[i],x2,y2[i]], /NoErase, CHARSIZE=plotchars,$
              thick=0.5
      index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
      if ctgd gt 0 then begin
         gdline = sort_linesall[index]
         gdlinelab = sort_linelaball[index]
         gdlinecol = sort_linecolall[index]
         gdlineyfrac = sort_lineyfracall[index]
         FOR M = 0,ctgd-1 DO BEGIN
            cgoplot,[gdline[M],gdline[M]],yran, color = gdlinecol[m], thick=1
            cgTEXT,gdline[m]-.1,yran[1]-(yran[0]+yran[1])*gdlineyfrac[m],$
                   gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
                   charthick=labchart,align=1
         endfor
      endif
         
   endfor

  CGPS_CLOSE
END
