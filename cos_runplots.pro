; 2016jul14  DSNR  created
pro cos_runplots

   table = '/Users/drupke/Dropbox/qsos/qsos.csv'
   fitdir = '/Users/drupke/Box Sync/cosquest/fits/'
   specdir = '/Users/drupke/Box Sync/cosquest/spectra/'
   plotdir = '/Users/drupke/Box Sync/cosquest/plots/spectra/'

   ifsf_fitdoublet,table,fitdir,'pg0844','NV'
;   ifsf_fitdoublet,table,fitdir,'pg0923','OVI'
   ifsf_fitdoublet,table,fitdir,'pg0953','OVI'
   ifsf_fitdoublet,table,fitdir,'pg1004','OVI'
   ifsf_fitdoublet,table,fitdir,'pg1116','OVI'
   ifsf_fitdoublet,table,fitdir,'pg1126','NV'
   ifsf_fitdoublet,table,fitdir,'pg1307','OVI'
   ifsf_fitdoublet,table,fitdir,'pg1351','NV'
   ifsf_fitdoublet,table,fitdir,'pg1411','NV'
   ifsf_fitdoublet,table,fitdir,'pg1440','NV'
   ifsf_fitdoublet,table,fitdir,'pg1613','NV'
   ifsf_fitdoublet,table,fitdir,'pg1613','OVI'
   ifsf_fitdoublet,table,fitdir,'pg1617','NV'
   ifsf_fitdoublet,table,fitdir,'pg1617','OVI'
   ifsf_fitdoublet,table,fitdir,'pg2130','NV'
   ifsf_fitdoublet,table,fitdir,'pg2214','NV'

   cos_uvlineplots,table,fitdir,plotdir,'pg0844',/NV
;   cos_uvlineplots,table,fitdir,plotdir,'pg0923',/OVI
   cos_uvlineplots,table,fitdir,plotdir,'pg0953',/OVI
   cos_uvlineplots,table,fitdir,plotdir,'pg1004',/OVI
   cos_uvlineplots,table,fitdir,plotdir,'pg1116',/OVI
   cos_uvlineplots,table,fitdir,plotdir,'pg1126',/NV
   cos_uvlineplots,table,fitdir,plotdir,'pg1307',/OVI
   cos_uvlineplots,table,fitdir,plotdir,'pg1351',/NV
   cos_uvlineplots,table,fitdir,plotdir,'pg1411',/NV
   cos_uvlineplots,table,fitdir,plotdir,'pg1440',/NV
   cos_uvlineplots,table,fitdir,plotdir,'pg1613',/NV,/OVI
   cos_uvlineplots,table,fitdir,plotdir,'pg1617',/NV,/OVI
   cos_uvlineplots,table,fitdir,plotdir,'pg2130',/NV
   cos_uvlineplots,table,fitdir,plotdir,'pg2214',/NV

   cos_spectraplots,table,specdir,plotdir,'izw1'
   cos_spectraplots,table,specdir,plotdir,'mrk231'
   cos_spectraplots,table,specdir,plotdir,'pg0007'
   cos_spectraplots,table,specdir,plotdir,'pg0026'
   cos_spectraplots,table,specdir,plotdir,'pg0157'
   cos_spectraplots,table,specdir,plotdir,'pg0804'
   cos_spectraplots,table,specdir,plotdir,'pg0838'
   cos_spectraplots,table,specdir,plotdir,'pg0844'
   cos_spectraplots,table,specdir,plotdir,'pg0923'
   cos_spectraplots,table,specdir,plotdir,'pg0953'
   cos_spectraplots,table,specdir,plotdir,'pg1004'
   cos_spectraplots,table,specdir,plotdir,'pg1116'
   cos_spectraplots,table,specdir,plotdir,'pg1126'
   cos_spectraplots,table,specdir,plotdir,'pg1229'
   cos_spectraplots,table,specdir,plotdir,'pg1302'
   cos_spectraplots,table,specdir,plotdir,'pg1307'
   cos_spectraplots,table,specdir,plotdir,'pg1309'
   cos_spectraplots,table,specdir,plotdir,'pg1351'
   cos_spectraplots,table,specdir,plotdir,'pg1411'
   cos_spectraplots,table,specdir,plotdir,'pg1435'
   cos_spectraplots,table,specdir,plotdir,'pg1440'
   cos_spectraplots,table,specdir,plotdir,'pg1613'
   cos_spectraplots,table,specdir,plotdir,'pg1617'
   cos_spectraplots,table,specdir,plotdir,'pg1626'
   cos_spectraplots,table,specdir,plotdir,'pg2130'
   cos_spectraplots,table,specdir,plotdir,'pg2214'
   cos_spectraplots,table,specdir,plotdir,'pg2349'

end
