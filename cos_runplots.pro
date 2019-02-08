; 2016jul14  DSNR  created
pro cos_runplots

   table = '/Users/drupke/Dropbox/qsos/qsos.csv'
   fitdir = '/Users/drupke/Box Sync/cosquest/fits/'
   specdir = '/Users/drupke/Box Sync/cosquest/spectra/'
   plotdir = '/Users/drupke/Box Sync/cosquest/plots/spectra/'

;   ifsf_fitdoublet,table,fitdir,'izw1','NV'
;   ifsf_fitdoublet,table,fitdir,'pg0804','NV'
;   ifsf_fitdoublet,table,fitdir,'pg0804','OVI'
;   ifsf_fitdoublet,table,fitdir,'pg0844','NV'
;;;   ifsf_fitdoublet,table,fitdir,'pg0923','OVI'
;   ifsf_fitdoublet,table,fitdir,'pg0953','OVI'
;   ifsf_fitdoublet,table,fitdir,'pg1001','OVI'
   ifsf_fitdoublet,table,fitdir,'pg1001','NV'
;   ifsf_fitdoublet,table,fitdir,'pg1004','OVI'
;   ifsf_fitdoublet,table,fitdir,'pg1116','OVI'
;   ifsf_fitdoublet,table,fitdir,'pg1126','NV'
;   ifsf_fitdoublet,table,fitdir,'pg1126','PV'
;   ifsf_fitdoublet,table,fitdir,'pg1307','OVI'
;   ifsf_fitdoublet,table,fitdir,'pg1351','NV'
;   ifsf_fitdoublet,table,fitdir,'pg1411','NV'
;   ifsf_fitdoublet,table,fitdir,'pg1411','PV'
;   ifsf_fitdoublet,table,fitdir,'pg1440','NV'
;   ifsf_fitdoublet,table,fitdir,'pg1448','NV'
;   ifsf_fitdoublet,table,fitdir,'pg1613','NV'
;   ifsf_fitdoublet,table,fitdir,'pg1613','OVI'
;   ifsf_fitdoublet,table,fitdir,'pg1617','NV'
;   ifsf_fitdoublet,table,fitdir,'pg1617','OVI'
;   ifsf_fitdoublet,table,fitdir,'pg2130','NV'
;   ifsf_fitdoublet,table,fitdir,'pg2214','NV'
;   ifsf_fitdoublet,table,fitdir,'pg2233','OVI'


;   cos_lineplots,table,fitdir,plotdir,'izw1',/NV,/lya,velran=[-3d3,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg0804',/NV,/OVI,/lya,velran=[-500,1.5d3]
;   cos_lineplots,table,fitdir,plotdir,'pg0844',/NV,/lya,velran=[-1.5d3,500]
;   cos_lineplots,table,fitdir,plotdir,'pg0953',/OVI,/lyb,/lyd,velran=[-1.5d3,500]
   cos_lineplots,table,fitdir,plotdir,'pg1001',/NV,/OVI,/lya,velran=[-1.3d4,2d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1004',/OVI,/lyb,/lyg,velran=[-1.2d4,0],$
;                 ignore=[1210,1220]
;   cos_lineplots,table,fitdir,plotdir,'pg1116',/OVI,/lya,/lyb,velran=[-6d3,0]
;   cos_lineplots,table,fitdir,plotdir,'pg1126',/NV,/PV,/lya,velran=[-5d3,1.5d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1307',/OVI,/lya,/lyb,velran=[-4.5d3,-2.5d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1351',/NV,/lya,velran=[-4d3,5d2]
;   cos_lineplots,table,fitdir,plotdir,'pg1411',/NV,/PV,/lya,velran=[-5d3,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1440',/NV,/lya,velran=[-3d3,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1448',/NV,/lya,velran=[-3d3,5d2]
;   cos_lineplots,table,fitdir,plotdir,'pg1613',/NV,/OVI,/lya,/lyb,velran=[-5d3,0]
;   cos_lineplots,table,fitdir,plotdir,'pg1617',/NV,/OVI,/lya,velran=[-3.5d3,0]
;   cos_lineplots,table,fitdir,plotdir,'pg2130',/NV,/lya,velran=[-3d3,1.5d3]
;   cos_lineplots,table,fitdir,plotdir,'pg2214',/NV, /lya, velran=[-4d3,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg2233',/OVI, /lyb, velran=[-2d3,1d3]

;   cos_spectraplots,table,specdir+'izw1.txt',$
;                    plotdir+'izw1_spec.eps','izw1'
;   cos_spectraplots,table,specdir+'mrk231.txt',$
;                    plotdir+'mrk231_spec.eps','mrk231',$
;                    ignore=[[1198,1202],[1301,1307],[1354,1360]]
;   cos_spectraplots,table,specdir+'pg0007.txt',$
;                    plotdir+'pg0007_spec.eps','pg0007'
;   cos_spectraplots,table,specdir+'pg0026.txt',$
;                    plotdir+'pg0026_spec.eps','pg0026'
;   cos_spectraplots,table,specdir+'pg0157.txt',$
;                    plotdir+'pg0157_spec.eps','pg0157'
;   cos_spectraplots,table,specdir+'pg0804.txt',$
;                    plotdir+'pg0804_spec.eps','pg0804'
;   cos_spectraplots,table,specdir+'pg0838.txt',$
;                    plotdir+'pg0838_spec.eps','pg0838',ignore=[1301,1307]
;   cos_spectraplots,table,specdir+'pg0844.txt',$
;                    plotdir+'pg0844_spec.eps','pg0844'
;   cos_spectraplots,table,specdir+'pg0923.txt',$
;                    plotdir+'pg0923_spec.eps','pg0923'
;   cos_spectraplots,table,specdir+'pg0953.txt',$
;                    plotdir+'pg0953_spec.eps','pg0953',ignore=[1301,1307]
;   cos_spectraplots,table,specdir+'pg1001.txt',$
;                    plotdir+'pg1001_spec.eps','pg1001',ignore=[1301,1307]
;   cos_spectraplots,table,specdir+'pg1004.txt',$
;                    plotdir+'pg1004_spec.eps','pg1004',ignore=[1448,1452]
;   cos_spectraplots,table,specdir+'pg1116.txt',$
;                    plotdir+'pg1116_spec.eps','pg1116',$
;                    ignore=[[1301,1307],[1430,1431]]
;   cos_spectraplots,table,specdir+'pg1126.txt',$
;                    plotdir+'pg1126_spec.eps','pg1126'
;   cos_spectraplots,table,specdir+'pg1126_cenwave1055.txt',$
;                    plotdir+'pg1126_cenwave1055_spec.eps','pg1126'
;   cos_spectraplots,table,specdir+'pg1126_withG160M.txt',$
;                    plotdir+'pg1126_withG160M_spec.eps','pg1126'
;   cos_spectraplots,table,specdir+'pg1202.txt',$
;                    plotdir+'pg1202_spec.eps','pg1202',ignore=[1301,1307]
;   cos_spectraplots,table,specdir+'pg1211.txt',$
;                    plotdir+'pg1211_spec.eps','pg1211'
;   cos_spectraplots,table,specdir+'pg1226.txt',$
;                    plotdir+'pg1226_spec.eps','pg1226'
;   cos_spectraplots,table,specdir+'pg1229.txt',$
;                    plotdir+'pg1229_spec.eps','pg1229'
;   cos_spectraplots,table,specdir+'pg1302.txt',$
;                    plotdir+'pg1302_spec.eps','pg1302',$
;                    ignore=[[1301,1307],[1457,1459]]
;   cos_spectraplots,table,specdir+'pg1307.txt',$
;                    plotdir+'pg1307_spec.eps','pg1307'
;   cos_spectraplots,table,specdir+'pg1309.txt',$
;                    plotdir+'pg1309_spec.eps','pg1309'
;   cos_spectraplots,table,specdir+'pg1351.txt',$
;                    plotdir+'pg1351_spec.eps','pg1351'
;   cos_spectraplots,table,specdir+'pg1411.txt',$
;                    plotdir+'pg1411_spec.eps','pg1411',ignore=[1301,1307]
;   cos_spectraplots,table,specdir+'pg1435.txt',$
;                    plotdir+'pg1435_spec.eps','pg1435'
;   cos_spectraplots,table,specdir+'pg1440.txt',$
;                    plotdir+'pg1440_spec.eps','pg1440'
;   cos_spectraplots,table,specdir+'pg1448.txt',$
;                    plotdir+'pg1448_spec.eps','pg1448'
;   cos_spectraplots,table,specdir+'pg1501.txt',$
;                    plotdir+'pg1501_spec.eps','pg1501'
;   cos_spectraplots,table,specdir+'pg1613.txt',$
;                    plotdir+'pg1613_spec.eps','pg1613',ignore=[1448,1450]
;   cos_spectraplots,table,specdir+'pg1617.txt',$
;                    plotdir+'pg1617_spec.eps','pg1617',ignore=[1450,1452]
;   cos_spectraplots,table,specdir+'pg1626.txt',$
;                    plotdir+'pg1626_spec.eps','pg1626',ignore=[1301,1307]
;   cos_spectraplots,table,specdir+'pg2130.txt',$
;                    plotdir+'pg2130_spec.eps','pg2130',ignore=[1458,1460]
;   cos_spectraplots,table,specdir+'pg2214.txt',$
;                    plotdir+'pg2214_spec.eps','pg2214',ignore=[1301,1307]
;   cos_spectraplots,table,specdir+'pg2233.txt',$
;                    plotdir+'pg2233_spec.eps','pg2233',ignore=[1301,1307]
;   cos_spectraplots,table,specdir+'pg2349.txt',$
;                    plotdir+'pg2349_spec.eps','pg2349'

end
