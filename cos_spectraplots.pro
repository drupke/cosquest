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
;    initfile: in, required, type=strarr
;      File that holds a list of galaxies as well as their redshifts.
;    directoryname: in, required, type=strarr
;      Location of the data files.
;    galaxyname: in, required, type=strarr
;      Galaxy being plotted, stored in a directory with the same name.
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
;
; :Copyright:
;    Copyright (C) 2016 Anthony To, David S. N. Rupke
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
PRO cos_spectraplots, table, specdir, plotdir, galaxyshortname


;List of emission/absorption lines and their corresponding wavelengths.

  LineWavelength = $
    [977.0201,$
    989.799,$
    1025.7223,$
    1031.912,$
    1037.613,$
    1122.5240,$
    1125.4477,$
    1128.0078,$
    1134.1653,$
    1134.4149,$
    1134.9803,$
    1143.2260,$
    1144.9379,$
    1152.8180,$
    1168.599,$
    1168.873,$
    1168.990,$
    1175.7217,$
    1183.030,$
    1184.544,$
    1190.416,$
    1193.2890,$
    1199.5496,$
    1200.2233,$
    1200.7098,$
    1206.500,$
    1215.67,$
    1238.8210,$
    1242.804,$
    1250.578,$
    1253.8051,$
    1259.518,$
    1260.4221,$
    1277.245,$
    1280.14,$
    1302.1685,$
    1304.3702,$
    1317.21,$
    1328.83,$
    1334.532,$
    1335.6627,$
    1335.7077,$
    1347.2396,$
    1355.5977,$
    1358.77,$
    1370.132,$
    1393.76,$
    1400.450,$
    1402.7729,$
    1548.202,$
    1550.774]
    
  LineLabel = $
    ['C III 977',$
    'N III 990',$
    'Ly$\beta$ 1026',$
    'O VI 1032',$
    'O VI 1038',$
    'Fe III 1123',$
    'Fe II 1025',$
    'P V 1128',$
    'N I 1134.2',$
    'N I 1134.4',$
    'N I 1134.9',$
    'Fe II 1143',$
    'Fe II 1145',$
    'P II 1153',$
    'N IV 1169',$
    'C IV 1169',$
    'C IV 1169',$
    'C III 1176',$
    'N III 1183',$
    'N III 1185',$
    'Si II 1190',$
    'Si II 1193',$
    'N I 1199.5',$
    'N I 1200.2',$
    'N I 1200.7',$
    'Si III 1207',$
    'Ly$\alpha$ 1216',$
    'N V 1239',$
    'N V 1243',$
    'S II 1251',$
    'S II 1254',$
    'S II 1260', $
    'Si II 1260',$
    'C I 1277',$
    'C I 1280',$
    'O I 1302',$
    'Si II 1304',$
    'Ni II 1317',$
    'C I 1329',$
    'C II 1335',$
    'C II* 1335.6',$
    'C II* 1335.7',$
    'C II 1347',$
    'O I 1356',$
    'Cu II 1359',$
    'Ni II 1370',$
    'Si IV 1394',$
    'Sn II 1400',$
    'Si IV 1403',$
    'C IV 1548',$
    'C IV 1551']
    

;Read wavelength and flux
  readcol, specdir+galaxyshortname+'.txt', wave, flux,$
           /silent

;Read galaxy full names and redshifts
  trows=[3,81]
  name = read_csvcol(table,'A',rows=trows,sep=',',type='string')
  galaxyshortnamelist = read_csvcol(table,'C',rows=trows,sep=',',type='string')
  z = read_csvcol(table,'D',rows=trows,sep=',',junk=bad)

  selectionparameter=WHERE(galaxyshortnamelist eq galaxyshortname)
  galaxyfullname=name[selectionparameter[0]]
  zsys=z[selectionparameter[0]]

;Set ASPECT RATIO for plots. Equivalent to (y-range)/(x-range) in data coords.
  aratio=(.11)/(.9)
  
;Shifts the absorption/emission waveelngths by the galaxy's redshift
  ShiftedLines= LineWavelength*(1+zsys)
  
;Rescaling flux values
  relativeflux=flux/1E-14
  
;Set the y-range of the big spectra
  yran=[0,Max(relativeflux)]

;Acquires the range of the wavelength values
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
  cgPS_OPEN,plotdir+galaxyshortname+'fullspectrum.eps',$
            /Encapsulated,scale_factor=1,Charsize=.5,/NOMATCH,$
            xsize=xsize, ysize=ysize, xoffset=xoffset, yoffset=yoffset, /Inches
;  !Y.MINOR=
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
    Position = [x1,.86,x2,.97], CHARSIZE=plotchars,thick=1,/NoErase,$
    xticklen=0.05,yticklen=0.01
    
;Plot legend
  AL_LEGEND,['Intrinsic','Galactic'], $
    Color=['Red','Blue'], charsize=1, charthick=2, $
    Linestyle=[0,0], $
    Position = [Min(wave)+.04*(Max(wave)-Min(wave)),$
                .9*(Max(relativeflux)-Min(relativeflux))], $
    bthick=.6,clear=1,linsize=0.25


;Plots 6 zoomed-in regions, one after another, along with absorption/emission lines and labels

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, $
    xran=[Min(wave),Min(wave)+(1d/6)*baserange], $
    axiscolor='Black',color='Black',yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(0d/6)*baserange): $
    value_locate(wave,Min(wave)+(1d/6)*baserange)])^2))], $
    Position = [x1,.695,x2,.805], /NoErase, CHARSIZE=plotchars,thick=1
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave),Min(wave)+(1d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(0d/6)*baserange): $
        value_locate(wave,Min(wave)+(1d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars, $
        charthick=labchart
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(0d/6)*baserange): $
        value_locate(wave,Min(wave)+(1d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars, $
        charthick=labchart
    ENDFOR
  ENDFOR

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, $
    xran=[Min(wave)+(1d/6)*baserange,Min(wave)+(2d/6)*baserange], $
    axiscolor='Black',color='Black',yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(1d/6)*baserange): $
    value_locate(wave,Min(wave)+(2d/6)*baserange)])^2))], $
    Position = [x1,.56,x2,.67], /NoErase, CHARSIZE=plotchars,$
    thick=1
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(1d/6)*baserange,Min(wave)+(2d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(1d/6)*baserange): $
        value_locate(wave,Min(wave)+(2d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars, $
        charthick=labchart
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(1d/6)*baserange): $
        value_locate(wave,Min(wave)+(2d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars, $
        charthick=labchart
    ENDFOR
  ENDFOR

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, xran=[Min(wave)+(2d/6)*baserange,Min(wave)+(3d/6)*baserange], $
    axiscolor='Black',color='Black',yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(2d/6)*baserange): $
    value_locate(wave,Min(wave)+(3d/6)*baserange)])^2))], $
    Position = [x1,.425,x2,.535], /NoErase, CHARSIZE=plotchars,$
    thick=1
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(2d/6)*baserange,Min(wave)+(3d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(2d/6)*baserange): $
        value_locate(wave,Min(wave)+(3d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars, $
        charthick=labchart
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(2d/6)*baserange): $
        value_locate(wave,Min(wave)+(3d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars, $
        charthick=labchart
    ENDFOR
  ENDFOR

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, $
    xran=[Min(wave)+(3d/6)*baserange,Min(wave)+(4d/6)*baserange], $
    axiscolor='Black',color='Black', yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(3d/6)*baserange): $
    value_locate(wave,Min(wave)+(4d/6)*baserange)])^2))], $
    Position = [x1,.29,x2,.40], /NoErase, CHARSIZE=plotchars,$
    thick=1
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(3d/6)*baserange,Min(wave)+(4d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(3d/6)*baserange): $
        value_locate(wave,Min(wave)+(4d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars, $
        charthick=labchart
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(3d/6)*baserange): $
        value_locate(wave,Min(wave)+(4d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars, $
        charthick=labchart
    ENDFOR
  ENDFOR

  cgplot, wave, relativeflux, xstyle=1, ystyle=1, $
    xran=[Min(wave)+(4d/6)*baserange,Min(wave)+(5d/6)*baserange], $
    axiscolor='Black',color='Black', yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(4d/6)*baserange): $
    value_locate(wave,Min(wave)+(5d/6)*baserange)])^2))], $
    Position = [x1,.155,x2,.265], /NoErase, CHARSIZE=plotchars,$
    thick=1
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(4d/6)*baserange,Min(wave)+(5d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(4d/6)*baserange): $
        value_locate(wave,Min(wave)+(5d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars,$
        charthick=labchart
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(4d/6)*baserange): $
        value_locate(wave,Min(wave)+(5d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars,$
        charthick=labchart
    ENDFOR
  ENDFOR
  
  cgplot, wave, relativeflux, xstyle=1, ystyle=1, $
    xran=[Min(wave)+(5d/6)*baserange,Min(wave)+(6d/6)*baserange], $
    axiscolor='Black',color='Black', yran=[0, $
    3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(5d/6)*baserange): $
    value_locate(wave,Min(wave)+(6d/6)*baserange)])^2))], $
    Position = [x1,.02,x2,.13], /NoErase, CHARSIZE=plotchars,$
    thick=1
  FOR M = 0, N_ELEMENTS(LineWavelength)-1 DO BEGIN
    cgoplot, [LineWavelength[M],LineWavelength[M]], 2*yran, color = 'Blue', thick=1
    cgoplot, [ShiftedLines[M],ShiftedLines[M]], 2*yran, color = 'Red', thick=1
    xran=[Min(wave)+(5d/6)*baserange,Min(wave)+(6d/6)*baserange]
    index=Where(LineWavelength lt xran[1] AND LineWavelength gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,LineWavelength[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(5d/6)*baserange): $
        value_locate(wave,Min(wave)+(6d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars,$
        charthick=labchart
    ENDFOR
    index=Where(ShiftedLines lt xran[1] AND ShiftedLines gt xran[0])
    FOR I = Min(index), Max(index) DO BEGIN
      cgTEXT,ShiftedLines[I]-.1, $
        .5*3*Sqrt(Mean((relativeflux[value_locate(wave,Min(wave)+(5d/6)*baserange): $
        value_locate(wave,Min(wave)+(6d/6)*baserange)])^2)),  $
        LineLabel[I], $
        /Data, $
        ORIENTATION = 90, $
        CHARSIZE = labchars,$
        charthick=labchart
    ENDFOR
  ENDFOR
  CGPS_CLOSE
END
