; docformat = 'rst'
;
;+
;
; Produces plots showing the normalized continuum with absorption profiles
; as well as the continuum fit.
;
; :Categories:
;    IFSFIT
;
; :Returns:
;    Postscript file with plots.
;
; :Params:
;    table: in, required, type=string
;      File that holds a list of galaxies and redshifts.
;    fitdir: in, required, type=string
;      Location of the data files.
;    galaxyshortname: in, required, type=string
;      Galaxy being plotted, stored in a fitdirectory with the same name.
;    plotdir: in, required, type=string
;      Location of the output plots.
;
; :Keywords:
;    ignorewave: in, optional, type=dblarr(N,2)
;      List of N regions to ignore in determining vertical plot range. Each
;      region has a lower and upper wavelength limit.
;    (AT LEAST ONE KEYWORD IS REQUIRED)
;    NV: in, optional, type=byte
;      Set if galaxy has been processed and NV doublets were found.
;    OVI: in, optional, type=byte
;      Set if galaxy has been processed and OVI doublets were found.
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
;      2016jul14, DSNR, cosmetic and input changes
;      2016sep07, DSNR, changed line label logic; moved line list to
;                       IFSF_LINELIST; added logic to prevent label collisions
;      2020aug05, DSNR, bug fix
;
; :Copyright:
;    Copyright (C) 2016--2020 Anthony To, David S. N. Rupke
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
PRO cos_lineplots, table, fitdir, plotdir, galaxyshortname, NV=NV, OVI=OVI, $
                   PV=PV,$
                   lya=lya, lyb=lyb, lyg=lyg, lyd=lyd, velran=velran,$
                   ignorewave=ignorewave, $
                   smooth=smooth,velnints=velnints

   c_kms = 299792.458d

;  List of emission/absorption lines and their corresponding wavelengths.
   linelab = 1b
   lines = ifsf_linelist(!NULL,linelab=linelab,/all,/vacuum)
;  Get lists from hashes
   keys = lines.keys()
   nlines = n_elements(keys)
   LineWavelength = dblarr(nlines)
   LineLabel = strarr(nlines)
   for i=0,nlines-1 do begin
      LineWavelength[i] = lines[keys[i]]
      LineLabel[i] = linelab[keys[i]]
   endfor

;  Read galaxy full names and redshifts
   trows=[3,85]
   name = read_csvcol(table,'A',rows=trows,sep=',',type='string')
   galaxyshortnamelist = read_csvcol(table,'C',rows=trows,sep=',',type='string')
   z = read_csvcol(table,'D',rows=trows,sep=',',junk=bad)

   selectionparameter=WHERE(galaxyshortnamelist eq galaxyshortname)
   galaxyfullname=name[selectionparameter[0]]
   zsys=z[selectionparameter[0]]

;  Shifts the emission/absorption wavelengths by the galaxy's redshift
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

;  Sets window sizes
   xpan_in = 3.5d
   ypan_in = 2.5d
   xmar_in = 0.6d
   ymar_in = 0.7d
   nline = keyword_set(NV) + keyword_set(OVI) + keyword_set(PV)
   if nline eq 3 then begin
      xsize=(xpan_in+xmar_in)*3d
      ysize=(ypan_in+ymar_in)*2d
      xmar_frac = xmar_in/xsize
      xpan_frac = xpan_in/xsize
      ymar_frac = ymar_in/ysize
      ypan_frac = ypan_in/ysize
;     Positions: top left, top right, bottom left, bottom right
      plotpos=dblarr(4,6)
      plotpos[*,0]=[xmar_frac,ymar_frac+ypan_frac,$
                    xmar_frac+xpan_frac,ymar_frac+2d*ypan_frac]
      plotpos[*,1]=[xmar_frac*2d + xpan_frac,ymar_frac+ypan_frac,$
                    xmar_frac*2d + xpan_frac*2d,ymar_frac+2d*ypan_frac]
      plotpos[*,2]=[xmar_frac*3d + xpan_frac*2d,ymar_frac+ypan_frac,$
                    1d,ymar_frac+2d*ypan_frac]
      plotpos[*,3]=[xmar_frac,ymar_frac,$
                    xmar_frac+xpan_frac,ymar_frac+ypan_frac]
      plotpos[*,4]=[xmar_frac*2d + xpan_frac,ymar_frac,$
                    xmar_frac*2d + xpan_frac*2d,ymar_frac+ypan_frac]
      plotpos[*,5]=[xmar_frac*3d + xpan_frac*2d,ymar_frac,$
                    1d,ymar_frac+ypan_frac]
   endif else if nline eq 2 then begin
      xsize=(xpan_in+xmar_in)*2d
      ysize=(ypan_in+ymar_in)*2d
      xmar_frac = xmar_in/xsize
      xpan_frac = xpan_in/xsize
      ymar_frac = ymar_in/ysize
      ypan_frac = ypan_in/ysize
;     Positions: top left, top right, bottom left, bottom right
      plotpos=dblarr(4,4)
      plotpos[*,0]=[xmar_frac,ymar_frac+ypan_frac,$
                    xmar_frac+xpan_frac,ymar_frac+2d*ypan_frac]
      plotpos[*,1]=[xmar_frac*2d + xpan_frac,ymar_frac+ypan_frac,$
                    1d,ymar_frac+2d*ypan_frac]
      plotpos[*,2]=[xmar_frac,ymar_frac,$
                    xmar_frac+xpan_frac,ymar_frac+ypan_frac]
      plotpos[*,3]=[xmar_frac*2d + xpan_frac,ymar_frac,$
                    1d,ymar_frac+ypan_frac]
   endif else begin
      xsize=xpan_in+xmar_in
      ysize=(ypan_in+ymar_in)*2d
      xmar_frac = xmar_in/xsize
      xpan_frac = xpan_in/xsize
      ymar_frac = ymar_in/ysize
      ypan_frac = ypan_in/ysize
;     Positions: top, bottom
      plotpos=dblarr(4,2)
      plotpos[*,0]=[xmar_frac,ymar_frac+ypan_frac,$
                    xmar_frac+xpan_frac,ymar_frac+2d*ypan_frac]
      plotpos[*,1]=[xmar_frac,ymar_frac,$
                    xmar_frac+xpan_frac,ymar_frac+ypan_frac]
   endelse

;  Legend position and sizes
   legpos = [xmar_frac*0.1d,1d - ymar_frac*0.1d]
   legcharsize = 0.6d
;  Title position and sizes
   if nline eq 1 then tit_xpos = xmar_frac+xpan_frac/2d $
   else if nline eq 2 then tit_xpos = xmar_frac*1.5d +xpan_frac $
   else tit_xpos = xmar_frac*2.5d +xpan_frac*2d
   tit_ypos = 1d - ymar_frac*0.5d
   titcharsize = 1.5d
;  Line label character size
   labchars=1d
  
;  Opens a Postscript file to draw plots on.

   CGPS_OPEN,plotdir+galaxyshortname+'_fit.eps',$
             /ENCAPSULATED, /NOMATCH,default_thick=2,charsize=1.5, xsize=xsize, $
             ysize=ysize, /Inches

   cgText,tit_xpos,tit_ypos,galaxyfullname +', '+ 'z='+String(zsys,format='(D0.3)'), $
          alignment=.5, Charsize = titcharsize
 
;  Track column
   icol = 0
 
;  PV

;  Read params
   IF keyword_set(PV) THEN BEGIN
      
      icol++
      
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'PV_fit_data.txt', $
               wave, modflux, continuum, flux, normalizedflux,format='(D,D,D,D,D)',$
               /silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'PV_fit_dataparam.txt', $
              xran_1,xran_2, format='(D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'PV_fit_datamodabs.txt', $
              nuvabs,elementsize,NUMLINE=1, format='(D,D)',/silent

;     Initializing and fixing variables  
      nuvabs=nuvabs[0]
      unity=make_array(elementsize,Value=1)
      elementsize=elementsize[0]
      moduvabs=Make_array(nuvabs,elementsize)

;     Rescaling flux values  
      continuum=continuum/1E-14
      flux = flux/1E-14
      moduvabs=moduvabs/1E-14

;     Defining plot ranges  
      xran=[xran_1[0],xran_2[0]]
      plotregion=where(wave ge xran_1[0] AND wave le xran_2[0])
;     Get rid of other regions
      if keyword_set(ignorewave) then begin
         signore = size(ignorewave)
         if signore[0] eq 1 then nignore=1 else nignore=signore[2]
         for i=0,nignore-1 do begin
            igd_tmp = where(wave le ignorewave[0,i] OR wave ge ignorewave[1,i])
            if i eq 0 then iignore_not = igd_tmp else $
               iignore_not = cgsetintersection(iignore_not,igd_tmp)
         endfor
         igdwave = cgsetintersection(iignore_not,plotregion)
      endif else igdwave = plotregion
  
;     Y-range for the normalized spectra plot
;      yrannormalized=[0,1.05*Max(normalizedflux[igdwave])]
      yrannormalized=[0,1.5d]
  
;     Y-range for the continuum plot
      yran=[0,1.05*Max(flux[igdwave])]
  
;     Normalized plot and legend
      normalizedfluxuse = normalizedflux
      if keyword_set(smooth) then $
         if smooth.haskey('PV') then $
            normalizedfluxuse = smooth(normalizedflux,smooth['PV'])
      cgplot,wave,normalizedfluxuse,xran=xran,yran=yrannormalized,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             xtickformat="(A1)",ytit='Normalized F!I$\lambda$',$
             position=plotpos[*,0], CHARSIZE=1,thick=1, Title= 'P V',/noerase
      AL_LEGEND,['Intrinsic','Galactic','Continuum','Components','Continuum w/Absorption'], $
                Color=['Red','Blue','Orange','Sky Blue','Purple'], charsize=legcharsize, $
                Linestyle=[0,0,0,0,0], Position =legpos,/norm,box=0
    
;     Plots absorption features, unity, and outlines spectra
      FOR M=0, nuvabs-1 DO BEGIN
         readcol, fitdir+'/'+galaxyshortname+'/'+galaxyshortname+$
                 'PV_fit_datamodabs.txt',moduvabs,skipline=1+elementsize*M, $
                 NUMLINE=elementsize, format='(D)',/silent
         cgoplot,wave,moduvabs,color='Sky Blue',thick=2
      ENDFOR
;      cgoplot,wave,unity,color='Orange',thick=2
      cgoplot,wave,modflux, color = 'Purple', thick = 4
  
;     Plots absorption and emission labels
      index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
      if ctgd gt 0 then begin
         gdline = sort_linesall[index]
         gdlinelab = sort_linelaball[index]
         gdlinecol = sort_linecolall[index]
         gdlineyfrac = sort_lineyfracall[index]
         FOR M = 0,ctgd-1 DO BEGIN
            cgoplot,[gdline[M],gdline[M]],yrannormalized,$
                    color = gdlinecol[m], thick=2
            cgTEXT,gdline[m]-.1,yrannormalized[1]-$
                   (yrannormalized[0]+yrannormalized[1])*gdlineyfrac[m],$
                   gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
                   charthick=labchart,align=1
         endfor
      endif

;     Plots continuum 
      if nline eq 3 then postmp = plotpos[*,3] $
      else if nline eq 2 then postmp = plotpos[*,2] $
      else postmp = plotpos[*,1]
      cgplot,wave,flux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             xtit='Observed Wavelength ($\Angstrom$)',$
             ytit='F!I$\lambda$!N/10!E-14!N (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)',$
             position=postmp,/NoErase, CHARSIZE=1,thick=1
      cgoplot,wave, continuum, color = 'Orange', thick=4
  
;     Plots absorption and emission labels
      index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
      if ctgd gt 0 then begin
         gdline = sort_linesall[index]
         gdlinelab = sort_linelaball[index]
         gdlinecol = sort_linecolall[index]
         gdlineyfrac = sort_lineyfracall[index]
         FOR M = 0,ctgd-1 DO BEGIN
            cgoplot,[gdline[M],gdline[M]],yran, color = gdlinecol[m], thick=2
            cgTEXT,gdline[m]-.1,yran[1]-$
                   (yran[0]+yran[1])*gdlineyfrac[m],$
                   gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
                   charthick=labchart,align=1
         endfor
      endif
  
   ENDIF 

 
;  OVI

;  Read params
   IF keyword_set(OVI) THEN BEGIN
      
      icol++

      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'OVI_fit_data.txt', $
               wave, modflux, continuum, flux, normalizedflux,format='(D,D,D,D,D)'
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'OVI_fit_dataparam.txt', $
              xran_1,xran_2, format='(D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'OVI_fit_datamodabs.txt', $
              nuvabs,elementsize,NUMLINE=1, format='(D,D)',/silent

;     Initializing and fixing variables  
      nuvabs=nuvabs[0]
      unity=make_array(elementsize,Value=1)
      elementsize=elementsize[0]
      moduvabs=Make_array(nuvabs,elementsize)

;     Rescaling flux values  
      continuum=continuum/1E-14
      flux = flux/1E-14
      moduvabs=moduvabs/1E-14

;     Defining plot ranges  
      xran=[xran_1[0],xran_2[0]]
      plotregion=where(wave ge xran_1[0] AND wave le xran_2[0])
;     Get rid of other regions
      if keyword_set(ignorewave) then begin
         signore = size(ignorewave)
         if signore[0] eq 1 then nignore=1 else nignore=signore[2]
         for i=0,nignore-1 do begin
            igd_tmp = where(wave le ignorewave[0,i] OR wave ge ignorewave[1,i])
            if i eq 0 then iignore_not = igd_tmp else $
               iignore_not = cgsetintersection(iignore_not,igd_tmp)
         endfor
         igdwave = cgsetintersection(iignore_not,plotregion)
      endif else igdwave = plotregion
  
;     Y-range for the normalized spectra plot
;      yrannormalized=[0,1.05*Max(normalizedflux[igdwave])]
      yrannormalized=[0,1.5d]

;     Y-range for the continuum plot
      yran=[0,1.05*Max(flux[igdwave])]
  
;     Normalized plot and legend
      postmp = plotpos[*,0]
      if nline eq 3 then postmp = plotpos[*,1] $
      else if nline eq 2 AND icol eq 2 then postmp = plotpos[*,1]
      normalizedfluxuse = normalizedflux
      if keyword_set(smooth) then $
         if smooth.haskey('OVI') then $
            normalizedfluxuse = smooth(normalizedflux,smooth['OVI'])
      cgplot,wave,normalizedfluxuse,xran=xran,yran=yrannormalized,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             xtickformat="(A1)",ytit='Normalized F!I$\lambda$',$
             position=postmp, CHARSIZE=1,thick=1, Title= 'O VI',/noerase
      AL_LEGEND,['Intrinsic','Galactic','Continuum','Components','Continuum w/Absorption'], $
                Color=['Red','Blue','Orange','Sky Blue','Purple'], charsize=legcharsize, $
                Linestyle=[0,0,0,0,0], Position =legpos,/norm,box=0
    
;     Plots absorption features, unity, and outlines spectra
      FOR M=0, nuvabs-1 DO BEGIN
         readcol, fitdir+'/'+galaxyshortname+'/'+galaxyshortname+$
                 'OVI_fit_datamodabs.txt',moduvabs,skipline=1+elementsize*M, $
                 NUMLINE=elementsize, format='(D)',/silent
         cgoplot,wave,moduvabs,color='Sky Blue',thick=2
      ENDFOR
;      cgoplot,wave,unity,color='Orange',thick=2
      cgoplot,wave,modflux, color = 'Purple', thick = 4
  
;     Plots absorption and emission labels
      index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
      if ctgd gt 0 then begin
         gdline = sort_linesall[index]
         gdlinelab = sort_linelaball[index]
         gdlinecol = sort_linecolall[index]
         gdlineyfrac = sort_lineyfracall[index]
         FOR M = 0,ctgd-1 DO BEGIN
            cgoplot,[gdline[M],gdline[M]],yrannormalized,$
                    color = gdlinecol[m], thick=2
            cgTEXT,gdline[m]-.1,yrannormalized[1]-$
                   (yrannormalized[0]+yrannormalized[1])*gdlineyfrac[m],$
                   gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
                   charthick=labchart,align=1
         endfor
      endif

  
;     Plots continuum 
      if nline eq 3 then postmp = plotpos[*,4] $
      else if nline eq 2 then begin
         if icol eq 2 AND icol eq 2 then postmp = plotpos[*,3] else postmp = plotpos[*,2]
      endif else postmp = plotpos[*,1]
      cgplot,wave,flux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             xtit='Observed Wavelength ($\Angstrom$)',$
             ytit='F!I$\lambda$!N/10!E-14!N (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)',$
             position=postmp,/NoErase, CHARSIZE=1,thick=1
      cgoplot,wave, continuum, color = 'Orange', thick=4
  
;     Plots absorption and emission labels
      index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
      if ctgd gt 0 then begin
         gdline = sort_linesall[index]
         gdlinelab = sort_linelaball[index]
         gdlinecol = sort_linecolall[index]
         gdlineyfrac = sort_lineyfracall[index]
         FOR M = 0,ctgd-1 DO BEGIN
            cgoplot,[gdline[M],gdline[M]],yran, color = gdlinecol[m], thick=2
            cgTEXT,gdline[m]-.1,yran[1]-$
                   (yran[0]+yran[1])*gdlineyfrac[m],$
                   gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
                   charthick=labchart,align=1
         endfor
      endif
  
   ENDIF 


;  NV

;  Read params
   IF keyword_set(NV) THEN BEGIN
      
      icol++
      
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'NV_fit_data.txt', $
              wave, modflux, continuum, flux, normalizedflux, format='(D,D,D,D,D)',$
              /silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'NV_fit_dataparam.txt', $
              xran_1,xran_2, format='(D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'NV_fit_datamodabs.txt', $
              nuvabs,elementsize,NUMLINE=1, format='(D,D)',/silent

;     Initializing and fixing params  
      nuvabs=nuvabs[0]
      unity=make_array(elementsize,Value=1)
      elementsize=elementsize[0]
      moduvabs=Make_array(nuvabs,elementsize)

;     Rescaling flux values  
      continuum=continuum/1E-14
      flux = flux/1E-14
      moduvabs=moduvabs/1E-14

;     Defining plot ranges  
      xran=[xran_1,xran_2]
      plotregion=where(wave ge xran_1[0] AND wave le xran_2[0])
;     Get rid of other regions
      if keyword_set(ignorewave) then begin
         signore = size(ignorewave)
         if signore[0] eq 1 then nignore=1 else nignore=signore[2]
         for i=0,nignore-1 do begin
            igd_tmp = where(wave le ignorewave[0,i] OR wave ge ignorewave[1,i])
            if i eq 0 then iignore_not = igd_tmp else $
               iignore_not = cgsetintersection(iignore_not,igd_tmp)
         endfor
         igdwave = cgsetintersection(iignore_not,plotregion)
      endif else igdwave = plotregion
  
;     Y-range for the normalized spectra plot
;      yrannormalized=[0,1.05*Max(normalizedflux[igdwave])]
      yrannormalized=[0,1.5d]

;     Y-range for the continuum plot
      yran=[0,1.05*Max(flux[igdwave])]
  
;     Normalized plot and legend
      if nline eq 3 then postmp = plotpos[*,2] $
      else if nline eq 2 then postmp = plotpos[*,1] $
      else postmp = plotpos[*,0]
      normalizedfluxuse = normalizedflux
      if keyword_set(smooth) then $
         if smooth.haskey('NV') then $
            normalizedfluxuse = smooth(normalizedflux,smooth['NV'])
      cgplot,wave,normalizedfluxuse,xran=xran,yran=yrannormalized,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             ytit='Normalized F!I$\lambda$',xtickformat="(A1)",$
             position=postmp, CHARSIZE=1,thick=1, Title= 'N V',/noerase
      AL_LEGEND,['Intrinsic','Galactic','Continuum','Components','Continuum w/Absorption'], $
                Color=['Red','Blue','Orange','Sky Blue','Purple'], charsize=legcharsize, $
                Linestyle=[0,0,0,0,0], Position=legpos,/norm,box=0
    
;     Plots absorption features, unity, and outlines spectra
      FOR M=0, nuvabs-1 DO BEGIN
         readcol, fitdir+'/'+galaxyshortname+'/'+galaxyshortname+'NV_fit_datamodabs.txt',$
                  moduvabs,skipline=1+elementsize*M, NUMLINE=elementsize, $
                  format='(D)',/silent
         cgoplot,wave,moduvabs,color='Sky Blue',thick=2
      ENDFOR
;      cgoplot,wave,unity,color='Orange',thick=2
      cgoplot,wave, modflux, color = 'Purple', thick = 4
  
;     Plots absorption and emission labels
      index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
      if ctgd gt 0 then begin
         gdline = sort_linesall[index]
         gdlinelab = sort_linelaball[index]
         gdlinecol = sort_linecolall[index]
         gdlineyfrac = sort_lineyfracall[index]
         FOR M = 0,ctgd-1 DO BEGIN
            cgoplot,[gdline[M],gdline[M]],yrannormalized,$
                    color = gdlinecol[m], thick=2
            cgTEXT,gdline[m]-.1,yrannormalized[1]-$
                   (yrannormalized[0]+yrannormalized[1])*gdlineyfrac[m],$
                   gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
                   charthick=labchart,align=1
         endfor
      endif

;     Plots continuum  
      if nline eq 3 then postmp = plotpos[*,5] $
      else if nline eq 2 then postmp = plotpos[*,3] $
      else postmp = plotpos[*,1]
      cgplot,wave,flux,xran=xran,yran=yran,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             ytit='F!I$\lambda$ !N/10!E-14!N (ergs s$\up-1$ cm$\up-2$ $\Angstrom$$\up-1$)',$
             position=postmp, /NoErase, CHARSIZE=1,thick=1,$
             xtit='Observed Wavelength ($\Angstrom$)'
      cgoplot,wave, continuum, color = 'Orange', thick=4
  
;     Plots absorption and emission labels
      index=Where(sort_linesall lt xran[1] AND sort_linesall gt xran[0],ctgd)
      if ctgd gt 0 then begin
         gdline = sort_linesall[index]
         gdlinelab = sort_linelaball[index]
         gdlinecol = sort_linecolall[index]
         gdlineyfrac = sort_lineyfracall[index]
         FOR M = 0,ctgd-1 DO BEGIN
            cgoplot,[gdline[M],gdline[M]],yran, color = gdlinecol[m], thick=2
            cgTEXT,gdline[m]-.1,yran[1]-$
                   (yran[0]+yran[1])*gdlineyfrac[m],$
                   gdlinelab[m],/Data,ORIENT=90d,CHARSIZE = labchars,$
                   charthick=labchart,align=1
         endfor
      endif
  
   ENDIF

   CGPS_CLOSE


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;  Sets window sizes
   xpan_in = 3.5d
   ypan_in = 1.5d
   xmar_in = 0.6d
   ymar_in = 0.7d
   nline = 0
   lys=[]
   if keyword_set(PV) then nline+=2
   if keyword_set(OVI) then nline+=2
   if keyword_set(NV) then nline+=2
   if keyword_set(lya) then begin
      nline++
      lys=[lys,'Lyalpha']
   endif
   if keyword_set(lyb) then begin
      nline++
      lys=[lys,'Lybeta']
   endif
   if keyword_set(lyg) then begin
      nline++
      lys=[lys,'Lygamma']
   endif
   if keyword_set(lyd) then begin
      nline++
      lys=[lys,'Lydelta']
   endif
   xsize=xpan_in+1.5d*xmar_in
   ysize=ypan_in*double(nline) + 1.5d*ymar_in
   xmar_frac = xmar_in/xsize
   xpan_frac = xpan_in/xsize
   ymar_frac = ymar_in/ysize
   ypan_frac = ypan_in/ysize
;  Positions, from bottom to top
   plotpos=dblarr(4,nline)
   for i=0,nline-1 do $
      plotpos[*,i]=[xmar_frac,ymar_frac+double(i)*ypan_frac,$
                    xmar_frac+xpan_frac,ymar_frac+double(i+1)*ypan_frac]
;  Title position and sizes
   tit_xpos = xmar_frac+0.5d*xpan_frac
   tit_ypos = 1d - ymar_frac*0.35d
   titcharsize = 1.5d
  
;  Opens a Postscript file to draw plots on.

   CGPS_OPEN,plotdir+galaxyshortname+'_vel.eps',$
             /ENCAPSULATED, /NOMATCH,default_thick=2,charsize=1.5, xsize=xsize, $
             ysize=ysize, /Inches

   if ~ keyword_set(velran) then velran=[-1d4,2d3]
   if ~ keyword_set(velnints) then velnints=0

   cgText,tit_xpos,tit_ypos,galaxyfullname +', '+ 'z='+String(zsys,format='(D0.3)'), $
          alignment=0.5, Charsize = titcharsize,/norm

   iplot = 0
   IF keyword_set(PV) THEN BEGIN
      linetmp = 'PV'
      blue='PV1117'
      red='PV1128'
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
              '_fit_data.txt',wave, modflux, continuum, flux, normalizedflux,$
              format='(D,D,D,D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
              '_fit_datamodabs.txt', $
              nuvabs,elementsize,NUMLINE=1, format='(D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
              'par_best.txt',veluvabs,skip=3,format='(X,X,X,X,D)',/silent

;     Initializing and fixing variables  
      nuvabs=nuvabs[0]
      elementsize=elementsize[0]
      moduvabs=dblarr(nuvabs,elementsize)

;     blue line
      zdiff = wave/(lines[blue]*(1d + zsys)) - 1d
      vel = c_kms * ((zdiff+1d)^2d - 1d) / $
                    ((zdiff+1d)^2d + 1d)
;      plotreg=[value_locate(vel,velran[0]):value_locate(vel,velran[1])]
;      yran=[0,1.05d*Max(normalizedflux[plotreg])]
      yran = [0,1.1d]
      normalizedfluxuse = normalizedflux
      if keyword_set(smooth) then $
         if smooth.haskey('PV') then $
            normalizedfluxuse = smooth(normalizedflux,smooth['PV'])
      cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             ytit='Normalized F!I$\lambda$',xtit='Velocity (km/s)',$
             position=plotpos[*,iplot], CHARSIZE=1,thick=1,/noerase,$
             xticks=velnints
      FOR M=0, nuvabs-1 DO BEGIN
         readcol, fitdir+'/'+galaxyshortname+'/'+galaxyshortname+$
                 linetmp+'_fit_datamodabs.txt',moduvabstmp,skipline=1+elementsize*M, $
                 NUMLINE=elementsize, format='(D)',/silent
         moduvabs[m,*] = moduvabstmp
         cgoplot,vel,moduvabstmp,color='Sky Blue',thick=2
         cgoplot,[veluvabs[m],veluvabs[m]],yran,thick=2,/linesty,$
                 color='Orange'
      ENDFOR
      cgoplot,vel,modflux, color = 'Purple', thick = 4
      xlinelab=velran[0]+0.02d*(velran[1]-velran[0])
      ylinelab=yran[0]+0.05*(yran[1]-yran[0])
      cgtext,linelab[blue],xlinelab,ylinelab,chars=1d,align=0
      iplot++

;     red line
      zdiff = wave/(lines[red]*(1d + zsys)) - 1d
      vel = c_kms * ((zdiff+1d)^2d - 1d) / $
                        ((zdiff+1d)^2d + 1d)
      yran = [0,1.1d]
      cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             xtickformat="(A1)",position=plotpos[*,iplot],$
             CHARSIZE=1,thick=1,/noerase,$
             xticks=velnints
      FOR M=0, nuvabs-1 DO begin
         cgoplot,vel,moduvabs[m,*],color='Sky Blue',thick=2
         cgoplot,[veluvabs[m],veluvabs[m]],yran,thick=2,/linesty,color='Orange'
      ENDFOR
      cgoplot,vel,modflux, color = 'Purple', thick = 4
      cgtext,linelab[red],xlinelab,ylinelab,chars=1d,align=0
      iplot++
      
      PVnuvabs = nuvabs
      PVveluvabs = veluvabs
  
   ENDIF 
   
   IF keyword_set(OVI) THEN BEGIN
      linetmp = 'OVI'
      blue='OVI1031'
      red='OVI1037'
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
              '_fit_data.txt',wave, modflux, continuum, flux, normalizedflux,$
              format='(D,D,D,D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
              '_fit_datamodabs.txt', $
              nuvabs,elementsize,NUMLINE=1, format='(D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
              'par_best.txt',veluvabs,skip=3,format='(X,X,X,X,D)',/silent
      modwave_ovi = wave
      modflux_ovi = modflux

;     Initializing and fixing variables  
      nuvabs=nuvabs[0]
      elementsize=elementsize[0]
      moduvabs=dblarr(nuvabs,elementsize)

;     blue line
      zdiff = wave/(lines[blue]*(1d + zsys)) - 1d
      vel = c_kms * ((zdiff+1d)^2d - 1d) / $
                    ((zdiff+1d)^2d + 1d)
;      plotreg=[value_locate(vel,velran[0]):value_locate(vel,velran[1])]
;      yran=[0,1.05d*Max(normalizedflux[plotreg])]
      yran = [0,1.1d]
      normalizedfluxuse = normalizedflux
      if keyword_set(smooth) then $
         if smooth.haskey('OVI') then $
            normalizedfluxuse = smooth(normalizedflux,smooth['OVI'])
      if keyword_set(PV) then $
         cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
                backg='White',axiscolor='Black',color='Black',$
                xtickformat="(A1)",$
                position=plotpos[*,iplot], CHARSIZE=1,thick=1,/noerase,$
                xticks=velnints $
      else $
         cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
                backg='White',axiscolor='Black',color='Black',$
                ytit='Normalized F!I$\lambda$',xtit='Velocity (km/s)',$
                position=plotpos[*,iplot], CHARSIZE=1,thick=1,/noerase,$
                xticks=velnints
       FOR M=0, nuvabs-1 DO BEGIN
         readcol, fitdir+'/'+galaxyshortname+'/'+galaxyshortname+$
                 linetmp+'_fit_datamodabs.txt',moduvabstmp,skipline=1+elementsize*M, $
                 NUMLINE=elementsize, format='(D)',/silent
         moduvabs[m,*] = moduvabstmp
         cgoplot,vel,moduvabstmp,color='Sky Blue',thick=2
         cgoplot,[veluvabs[m],veluvabs[m]],yran,thick=2,/linesty,$
                 color='Blue'
      ENDFOR
      cgoplot,vel,modflux, color = 'Purple', thick = 4
      xlinelab=velran[0]+0.02d*(velran[1]-velran[0])
      ylinelab=yran[0]+0.05*(yran[1]-yran[0])
      cgtext,linelab[blue],xlinelab,ylinelab,chars=1d,align=0
      iplot++

;     red line
      zdiff = wave/(lines[red]*(1d + zsys)) - 1d
      vel = c_kms * ((zdiff+1d)^2d - 1d) / $
                        ((zdiff+1d)^2d + 1d)
      yran = [0,1.1d]
      cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
             backg='White',axiscolor='Black',color='Black',$
             xtickformat="(A1)",position=plotpos[*,iplot],$
             CHARSIZE=1,thick=1,/noerase,$
             xticks=velnints
      FOR M=0, nuvabs-1 DO begin
         cgoplot,vel,moduvabs[m,*],color='Sky Blue',thick=2
         cgoplot,[veluvabs[m],veluvabs[m]],yran,thick=2,/linesty,color='Blue'
      ENDFOR
      cgoplot,vel,modflux, color = 'Purple', thick = 4
      cgtext,linelab[red],xlinelab,ylinelab,chars=1d,align=0
      iplot++
      
      OVInuvabs = nuvabs
      OVIveluvabs = veluvabs
  
   ENDIF 

   IF keyword_set(NV) THEN BEGIN
      linetmp = 'NV'
      blue='NV1238'
      red='NV1242'
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
         '_fit_data.txt',wave, modflux, continuum, flux, normalizedflux,$
         format='(D,D,D,D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
         '_fit_datamodabs.txt', $
         nuvabs,elementsize,NUMLINE=1, format='(D,D)',/silent
      readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
         'par_best.txt',veluvabs,skip=3,format='(X,X,X,X,D)',/silent
      modwave_nv = wave
      modflux_nv = modflux

;     Initializing and fixing variables
      nuvabs=nuvabs[0]
      elementsize=elementsize[0]
      moduvabs=dblarr(nuvabs,elementsize)

;     blue line
      zdiff = wave/(lines[blue]*(1d + zsys)) - 1d
      vel = c_kms * ((zdiff+1d)^2d - 1d) / $
                    ((zdiff+1d)^2d + 1d)
      yran = [0,1.1d]
      normalizedfluxuse = normalizedflux
      if keyword_set(smooth) then $
         if smooth.haskey('NV') then $
            normalizedfluxuse = smooth(normalizedflux,smooth['NV'])
      if keyword_set(OVI) OR keyword_set(PV) then $
         cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
                backg='White',axiscolor='Black',color='Black',$
                xtickformat="(A1)",$
                position=plotpos[*,iplot], CHARSIZE=1,thick=1,/noerase,$
                xticks=velnints $
      else $
         cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
                backg='White',axiscolor='Black',color='Black',$
                ytit='Normalized F!I$\lambda$',xtit='Velocity (km/s)',$
                position=plotpos[*,iplot], CHARSIZE=1,thick=1,/noerase,$
                xticks=velnints
      FOR M=0, nuvabs-1 DO BEGIN
         readcol, fitdir+'/'+galaxyshortname+'/'+galaxyshortname+$
            linetmp+'_fit_datamodabs.txt',moduvabstmp,skipline=1+elementsize*M, $
            NUMLINE=elementsize, format='(D)',/silent
         moduvabs[m,*] = moduvabstmp
         cgoplot,vel,moduvabstmp,color='Sky Blue',thick=2
         cgoplot,[veluvabs[m],veluvabs[m]],yran,thick=2,/linesty,$
                 color='Red'
      ENDFOR
      cgoplot,vel,modflux, color = 'Purple', thick = 4
      xlinelab=velran[0]+0.02d*(velran[1]-velran[0])
      ylinelab=yran[0]+0.05*(yran[1]-yran[0])
      cgtext,linelab[blue],xlinelab,ylinelab,chars=1d,align=0
     iplot++

;     red line
      zdiff = wave/(lines[red]*(1d + zsys)) - 1d
      vel = c_kms * ((zdiff+1d)^2d - 1d) / $
         ((zdiff+1d)^2d + 1d)
      yran = [0,1.1d]
      cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
         backg='White',axiscolor='Black',color='Black',$
         xtickformat="(A1)",position=plotpos[*,iplot],$
         CHARSIZE=1,thick=1,/noerase,$
         xticks=velnints
      FOR M=0, nuvabs-1 DO begin
         cgoplot,vel,moduvabs[m,*],color='Sky Blue',thick=2
         cgoplot,[veluvabs[m],veluvabs[m]],yran,thick=2,/linesty,color='Red'
      ENDFOR
      cgoplot,vel,modflux, color = 'Purple', thick = 4
      cgtext,linelab[red],xlinelab,ylinelab,chars=1d,align=0
      iplot++

      NVnuvabs = nuvabs
      NVveluvabs = veluvabs

   ENDIF
   
   IF n_elements(lys) gt 0 then begin
      foreach linetmp,lys do begin
         divmod = 0b
         readcol,fitdir+'/'+galaxyshortname+'/'+galaxyshortname+linetmp+$
                 '_contfit.txt',wave, normalizedflux,$
                 format='(D,X,X,D,X)',/silent

         zdiff = wave/(lines[linetmp]*(1d + zsys)) - 1d
         vel = c_kms * ((zdiff+1d)^2d - 1d) / $
                       ((zdiff+1d)^2d + 1d)
         yran = [0d,1.1d]
         normalizedfluxuse = normalizedflux
         if keyword_set(smooth) then $
            if smooth.haskey(linetmp) then $
               normalizedfluxuse = smooth(normalizedflux,smooth[linetmp])
         fluxdivmod = normalizedfluxuse


         cgplot,vel,normalizedfluxuse,xran=velran,yran=yran,xstyle=1,ystyle=1,$
            backg='White',axiscolor='Black',color='Black',$
            xtickformat="(A1)",$
            position=plotpos[*,iplot], CHARSIZE=1,thick=1,/noerase,$
            xticks=velnints


         if keyword_set(PV) then $
            FOR M=0, PVnuvabs-1 DO $
               cgoplot,[PVveluvabs[m],PVveluvabs[m]],yran,thick=2,/linesty,$
                       color='Orange'
         if keyword_set(OVI) then begin
            FOR M=0, OVInuvabs-1 DO $
               cgoplot,[OVIveluvabs[m],OVIveluvabs[m]],yran,thick=2,/linesty,$
                       color='Blue'
             if linetmp eq 'Lybeta' then begin
               divmod = 1b
               tmpmod = fluxdivmod / fluxdivmod
               ialign_left = where(wave eq modwave_ovi[0],ctalign)
               if ctalign gt 0 then begin
                  if n_elements(modflux_ovi) ge n_elements(fluxdivmod) - ialign_left then begin
                     ialign_right = n_elements(fluxdivmod)-1
                     jalign_right = ialign_right - ialign_left
                  endif else begin
                     ialign_right = ialign_left + n_elements(modflux_ovi)-1
                     jalign_right = n_elements(modflux_ovi)-1
                  endelse
                  fluxdivmod[ialign_left:ialign_right] /= modflux_ovi[0:jalign_right]
;  Trying to figure out noisy spectrum divided by very small model values ...
;                  fluxdivmodtmp = fluxdivmod[ialign_left:ialign_right]
;                  modfluxdiv = modflux_ovi[0:jalign_right]
;                  inearzero = where(modfluxdiv le 0.05,ctnearzero)
;                  if ctnearzero gt 0 then begin
;                     fluxdivmodtmp[inearzero] = 
;                     fluxdivmod[ialign_left:ialign_right] = fluxdivmodtmp
;                  endif
               endif else begin
                  print,'COS_LINEPLOTS: Cannot align OVI models with Lybeta data.'
               endelse
            endif
         endif
         if keyword_set(NV) then begin
            FOR M=0, NVnuvabs-1 DO $
               cgoplot,[NVveluvabs[m],NVveluvabs[m]],yran,thick=2,/linesty,$
                       color='Red'
             if linetmp eq 'Lyalpha' then begin
               divmod = 1b
               tmpmod = fluxdivmod / fluxdivmod
               ialign_left = where(wave eq modwave_nv[0],ctalign)
               if ctalign gt 0 then begin
                  if n_elements(modflux_nv) ge n_elements(fluxdivmod) - ialign_left then begin
                     ialign_right = n_elements(fluxdivmod)-1
                     jalign_right = ialign_right - ialign_left
                  endif else begin
                     ialign_right = ialign_left + n_elements(modflux_nv)-1
                     jalign_right = n_elements(modflux_nv)-1
                  endelse
                  fluxdivmod[ialign_left:ialign_right] /= modflux_nv[0:jalign_right]
               endif else begin
                  print,'COS_LINEPLOTS: Cannot align NV models with Lybeta data.'
               endelse
            endif
         endif
         
         if divmod then begin
            diff = normalizedfluxuse - fluxdivmod
            idiff = where(diff ne 0d,ctdiff)
            if ctdiff gt 0 then $
               cgoplot,vel[idiff],fluxdivmod[idiff],color='Purple'
         endif

         xlinelab=velran[0]+0.02d*(velran[1]-velran[0])
         ylinelab=yran[0]+0.05*(yran[1]-yran[0])
         cgtext,linelab[linetmp],xlinelab,ylinelab,chars=1d,align=0
         iplot++
         cgoplot,velran,[1d,1d]
      endforeach
   ENDIF

   CGPS_CLOSE

END
