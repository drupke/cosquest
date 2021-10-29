; 2016jul14  DSNR  created
pro cos_runplots

   table = '/Users/drupke/Box Sync/qsos/qsos.csv'
   fitdir = '/Users/drupke/Box Sync/cosquest/fits/'
   specdir = '/Users/drupke/Box Sync/cosquest/spectra/'
   plotdir = '/Users/drupke/Box Sync/cosquest/plots/spectra/'
   tabdir = '/Users/drupke/Box Sync/cosquest/tables/'

;   cos_izw1ly
;   cos_pg0804ly
;   cos_pg0844ly
;   cos_pg0923ly
;   cos_pg0953ly
;   cos_pg1001ly
;   cos_pg1004ly
;   cos_pg1116ly
;   cos_pg1126ly
;   cos_pg1307ly
;   cos_pg1309ly
;   cos_pg1351ly
;   cos_pg1411ly
;   cos_pg1440ly
;   cos_pg1448ly
;   cos_pg1613ly
;   cos_pg1617ly
;   cos_pg2130ly
;   cos_pg2214ly
;   cos_pg2233ly
;
   nsplit=8
   nomc=1b
;
   ifsf_fitdoublet,fitdir,'izw1','NV','cos_fitdoublet_galinfo',$
                      argsgalinfo={dir:fitdir,galshort:'izw1',doublet:'NV'},$
                      nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg0804','NV','cos_fitdoublet_galinfo',$
                      argsgalinfo={dir:fitdir,galshort:'pg0804',doublet:'NV'},$
                      nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg0804','OVI','cos_fitdoublet_galinfo',$
                      argsgalinfo={dir:fitdir,galshort:'pg0804',doublet:'OVI'},$
                      nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg0844','NV','cos_fitdoublet_galinfo',$
                      argsgalinfo={dir:fitdir,galshort:'pg0844',doublet:'NV'},$
                      nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg0923','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg0923',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg0953','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg0953',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1001','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1001',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1001','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1001',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1004','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1004',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1116','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1116',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1126','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1126',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1126','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1126',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1126','PV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1126',doublet:'PV'},$
                   nsplit=nsplit,nomc=nomc
;   ifsf_fitdoublet,fitdir,'mrk231','NV','cos_fitdoublet_galinfo',$
;                   argsgalinfo={dir:fitdir,galshort:'mrk231',doublet:'NV'},$
;                   /init
   ifsf_fitdoublet,fitdir,'pg1307','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1307',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1309','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1309',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1351','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1351',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1411','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1411',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1411','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1411',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1411','PV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1411',doublet:'PV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1440','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1440',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1448','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1448',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1613','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1613',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1613','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1613',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1617','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1617',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1617','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1617',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg1617','PV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg1617',doublet:'PV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg2130','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg2130',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg2214','NV','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg2214',doublet:'NV'},$
                   nsplit=nsplit,nomc=nomc
   ifsf_fitdoublet,fitdir,'pg2233','OVI','cos_fitdoublet_galinfo',$
                   argsgalinfo={dir:fitdir,galshort:'pg2233',doublet:'OVI'},$
                   nsplit=nsplit,nomc=nomc
;
;
;   cos_lineplots,table,fitdir,plotdir,'izw1',/NV,/lya,velran=[-2d3,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg0804',/NV,/OVI,/lya,velran=[0,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg0844',/NV,/lya,velran=[-500,500]
;   cos_lineplots,table,fitdir,plotdir,'pg0923',/OVI,/lya,/lyb,velran=[-4d3,0],$
;                 ignore=[1213,1219]
;   cos_lineplots,table,fitdir,plotdir,'pg0953',/OVI,/lyb,/lyd,velran=[-1.5d3,500]
;   cos_lineplots,table,fitdir,plotdir,'pg1001',/NV,/OVI,/lya,/lyb,/lyg,$
;      velran=[-1d4,1d3],smooth=hash("Lygamma",5d)
;   cos_lineplots,table,fitdir,plotdir,'pg1004',/OVI,/lyb,/lyg,velran=[-1.2d4,0],$
;                 ignore=[1210,1220],velnints=4
;   cos_lineplots,table,fitdir,plotdir,'pg1116',/OVI,/lya,/lyb,velran=[-4d3,0]
;   cos_lineplots,table,fitdir,plotdir,'pg1126',/OVI,/NV,/PV,/lya,/lyb,$
;                 velran=[-5d3,1.5d3],smooth=hash("OVI",5d,"Lybeta",10d)
;   cos_lineplots,table,fitdir,plotdir,'pg1307',/OVI,/lya,/lyb,velran=[-4.5d3,-2.5d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1309',/OVI,/lya,/lyb,/lyg,$
;                 velran=[-1.9d3,0.5d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1351',/NV,/lya,velran=[-4d3,5d2]
;   cos_lineplots,table,fitdir,plotdir,'pg1411',/NV,/OVI,/PV,/lya,/lyb,velran=[-5d3,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1440',/NV,/lya,velran=[-3d3,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg1448',/NV,/lya,velran=[-750d,250d]
;   cos_lineplots,table,fitdir,plotdir,'pg1613',/NV,/OVI,/lya,/lyb,velran=[-4500,-2500]
;   cos_lineplots,table,fitdir,plotdir,'pg1617',/NV,/OVI,/lya,/lyb,$
;      velran=[-5d3,0]
;   cos_lineplots,table,fitdir,plotdir,'pg2130',/NV,/lya,velran=[-2000d,500d]
;   cos_lineplots,table,fitdir,plotdir,'pg2214',/NV, /lya, velran=[-4d3,1d3]
;   cos_lineplots,table,fitdir,plotdir,'pg2233',/OVI, /lyb, velran=[-500d,0d]
;
;   medsn=1b
;   flam=1b
;   gallist = list()
;   medsnlist = list()
;   flamlist = list()
;   cos_spectraplots,table,specdir+'pg0007.txt',$
;                    plotdir+'pg0007_spec.eps','pg0007',medsn=medsn,flam=flam
;   gallist.add,'pg0007'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg0026.txt',$
;                    plotdir+'pg0026_spec.eps','pg0026',medsn=medsn,flam=flam
;   gallist.add,'pg0026'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'izw1.txt',$
;                    plotdir+'izw1_spec.eps','izw1',medsn=medsn,flam=flam
;   gallist.add,'izw1'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg0157.txt',$
;                    plotdir+'pg0157_spec.eps','pg0157',medsn=medsn,flam=flam
;   gallist.add,'pg0157'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg0804.txt',$
;                    plotdir+'pg0804_spec.eps','pg0804',medsn=medsn,flam=flam
;   gallist.add,'pg0804'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg0838.txt',$
;                    plotdir+'pg0838_spec.eps','pg0838',ignore=[1301,1307],medsn=medsn,flam=flam
;   gallist.add,'pg0838'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg0844.txt',$
;                    plotdir+'pg0844_spec.eps','pg0844',medsn=medsn,flam=flam
;   gallist.add,'pg0844'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg0923.txt',$
;                    plotdir+'pg0923_spec.eps','pg0923',medsn=medsn,flam=flam
;   gallist.add,'pg0923'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg0953.txt',$
;                    plotdir+'pg0953_spec.eps','pg0953',ignore=[1301,1307],medsn=medsn,flam=flam
;   gallist.add,'pg0953'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1001.txt',$
;                    plotdir+'pg1001_spec.eps','pg1001',ignore=[1066,1100],medsn=medsn,flam=flam
;   gallist.add,'pg1001'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1004.txt',$
;                    plotdir+'pg1004_spec.eps','pg1004',ignore=[1448,1452],medsn=medsn,flam=flam
;   gallist.add,'pg1004'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1116.txt',$
;                    plotdir+'pg1116_spec.eps','pg1116',$
;                    ignore=[[1301,1307],[1430,1431]],medsn=medsn,flam=flam
;   gallist.add,'pg1116'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1126.txt',$
;                    plotdir+'pg1126_spec.eps','pg1126',medsn=medsn,flam=flam
;   gallist.add,'pg1126'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1126_cenwave1055.txt',$
;                    plotdir+'pg1126_cenwave1055_spec.eps','pg1126',ignore=[900,1090]
;;   cos_spectraplots,table,specdir+'pg1126_withG160M.txt',$
;;                    plotdir+'pg1126_withG160M_spec.eps','pg1126'
;;   cos_spectraplots,table,specdir+'pg1202.txt',$
;;                    plotdir+'pg1202_spec.eps','pg1202',ignore=[1301,1307],medsn=medsn,flam=flam
;;   medsnlist.add, medsn
;;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1211.txt',$
;                    plotdir+'pg1211_spec.eps','pg1211',medsn=medsn,flam=flam
;   gallist.add,'pg1211'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1226.txt',$
;                    plotdir+'pg1226_spec.eps','pg1226',medsn=medsn,flam=flam
;   gallist.add,'pg1226'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1229.txt',$
;                    plotdir+'pg1229_spec.eps','pg1229',medsn=medsn,flam=flam
;   gallist.add,'pg1229'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'mrk231.txt',$
;                    plotdir+'mrk231_spec.eps','mrk231',$
;                    ignore=[[1198,1202],[1301,1307],[1354,1360]],medsn=medsn,flam=flam
;   gallist.add,'mrk231'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1302.txt',$
;                    plotdir+'pg1302_spec.eps','pg1302',$
;                    ignore=[[1301,1307],[1457,1459]],medsn=medsn,flam=flam
;   gallist.add,'pg1302'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1307.txt',$
;                    plotdir+'pg1307_spec.eps','pg1307',medsn=medsn,flam=flam
;   gallist.add,'pg1307'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1309.txt',$
;                    plotdir+'pg1309_spec.eps','pg1309',medsn=medsn,flam=flam
;   gallist.add,'pg1309'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1351.txt',$
;                    plotdir+'pg1351_spec.eps','pg1351',medsn=medsn,flam=flam
;   gallist.add,'pg1351'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1411.txt',$
;                    plotdir+'pg1411_spec.eps','pg1411',ignore=[1301,1307],medsn=medsn,flam=flam
;   gallist.add,'pg1411'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1411_cenwave1096.txt',$
;      plotdir+'pg1411_cenwave1096_spec.eps','pg1411'
;   cos_spectraplots,table,specdir+'pg1435.txt',$
;                    plotdir+'pg1435_spec.eps','pg1435',medsn=medsn,flam=flam
;   gallist.add,'pg1435'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1440.txt',$
;                    plotdir+'pg1440_spec.eps','pg1440',medsn=medsn,flam=flam
;   gallist.add,'pg1440'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1448.txt',$
;                    plotdir+'pg1448_spec.eps','pg1448',medsn=medsn,flam=flam
;   gallist.add,'pg1448'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1501.txt',$
;                    plotdir+'pg1501_spec.eps','pg1501',medsn=medsn,flam=flam
;   gallist.add,'pg1501'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1613.txt',$
;                    plotdir+'pg1613_spec.eps','pg1613',ignore=[1448,1450],medsn=medsn,flam=flam
;   gallist.add,'pg1613'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1617.txt',$
;                    plotdir+'pg1617_spec.eps','pg1617',ignore=[1450,1452],medsn=medsn,flam=flam
;   gallist.add,'pg1617'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg1626.txt',$
;                    plotdir+'pg1626_spec.eps','pg1626',ignore=[1301,1307],medsn=medsn,flam=flam
;   gallist.add,'pg1626'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg2130.txt',$
;                    plotdir+'pg2130_spec.eps','pg2130',ignore=[1458,1460],medsn=medsn,flam=flam
;   gallist.add,'pg2130'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg2214.txt',$
;                    plotdir+'pg2214_spec.eps','pg2214',ignore=[1301,1307],medsn=medsn,flam=flam
;   gallist.add,'pg2214'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg2233.txt',$
;                    plotdir+'pg2233_spec.eps','pg2233',ignore=[1301,1307],medsn=medsn,flam=flam
;   gallist.add,'pg2233'
;   medsnlist.add, medsn
;   flamlist.add,flam
;   cos_spectraplots,table,specdir+'pg2349.txt',$
;                    plotdir+'pg2349_spec.eps','pg2349',medsn=medsn,flam=flam
;   gallist.add,'pg2349'
;   medsnlist.add, medsn
;   flamlist.add,flam
;
;   medsnarr = medsnlist.toarray()
;   print,'S/N over sample at 1290-1310 A: '
;   print,'Median: ',median(medsnarr),format='(A0,D0.2)'
;   print,'Mean: ',mean(medsnarr),format='(A0,D0.2)'
;   print,'Min: ',min(medsnarr),format='(A0,D0.2)'
;   print,'Max: ',max(medsnarr),format='(A0,D0.2)'
;   print,'Std.dev.: ',stddev(medsnarr),format='(A0,D0.2)'
;
;   openw,lun_tmp,tabdir+'tab_specvals.txt',/get_lun
;   galarr = gallist.toarray()
;   flamarr = flamlist.toarray()
;   printf,lun_tmp,'#From COSQUEST spectra',format='(A0)
;   printf,lun_tmp,'#Flux units: 1e-14 erg/s/cm^2/A',format='(A0)
;   printf,lun_tmp,'#S/N is median over 20A band',format='(A0)
;   printf,lun_tmp,'#Gal','flam1125','snr1300',format='(A-8,A10,A10)
;;   printf,lun_tmp,'#Gal','flam1125','ferr1125','snr1300',format='(A-8,3A10)
;   for i=0,n_elements(galarr)-1 do $
;;      printf,lun_tmp,galarr[i],flamarr[i,0]/1d-14,flamarr[i,1]/1d-14,medsnarr[i],$
;;         format='(A-8,3D10.4)'
;      printf,lun_tmp,galarr[i],flamarr[i]/1d-14,medsnarr[i],$
;         format='(A-8,3D10.4)'
;   free_lun,lun_tmp


end
