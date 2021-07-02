; docformat = 'rst'
;
;+
;
; Initialize parameters for fitting. Specific to GMOS instrument.
;
; :Categories:
;    IFSFIT/INIT
;
; :Returns:
;    PARINFO structure for input into MPFIT.
;
; :Params:
;    linelist: in, required, type=hash(lines)
;      Emission line rest frame wavelengths.
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
;      2019nov18, DSNR, moved from IFSF_FITDOUBLET to stand-alone routine.
;    
; :Copyright:
;    Copyright (C) 2019 David S. N. Rupke
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
function cos_fitdoublet_galinfo,redshift,dir=dir,galshort=galshort,doublet=doublet

   table = '/Users/drupke/Box Sync/qsos/qsos.csv'
   trows=[3,85]

   galshortlist = read_csvcol(table,'C',rows=trows,sep=',',type='string')
   zlist = read_csvcol(table,'D',rows=trows,sep=',',junk=bad)
   
   selectionparameter=WHERE(galshortlist eq galshort)
   galshort=galshortlist[selectionparameter[0]]
   redshift=zlist[selectionparameter]
   readcol, dir+galshort+'/'+galshort+doublet+'par_init.txt', $
            profileshifts, profilesig, coveringfactor, opticaldepth, $
            FORMAT='(A,D,D,D,D)',/silent
   initproc = 'cos_'+galshort+doublet
   initstr = call_function(initproc,dir, galshort, redshift, $
                           profileshifts, profilesig, coveringfactor, $
                           opticaldepth)

   return,initstr

end
