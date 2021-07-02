pro cos_coadd

   datdir = '/Users/drupke/hst/cos/archivedownloads/'
   specdir = '/Users/drupke/Box Sync/cosquest/spectra/'
   fmt = '(D15.6,E15.6,E15.6)'

   ;PG0007:
   coadd_x1d,wave,flux,err,path=datdir+'pg0007/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg0007.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG0026:
   coadd_x1d,wave,flux,err,path=datdir+'pg0026/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg0026.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;I Zw 1:
   coadd_x1d,wave,flux,err,path=datdir+'izw1/',chan=1,plot=1,bin=3
   openw,lun,specdir+'izw1.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun

   ;PG0157:
   coadd_x1d,wave,flux,err,path=datdir+'pg0157/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg0157.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG0804:
   coadd_x1d,wave,flux,err,path=datdir+'pg0804/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg0804.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG0838:
   coadd_x1d,wave,flux,err,path=datdir+'pg0838/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg0838.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG0844:
   coadd_x1d,wave,flux,err,path=datdir+'pg0844/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg0844.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG0923:
   coadd_x1d,wave,flux,err,path=datdir+'pg0923/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg0923.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG0953:
   coadd_x1d,wave,flux,err,path=datdir+'pg0953/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg0953.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1001: Added Cooke CENWAVE=1222 exposure
   coadd_x1d,wave,flux,err,path=datdir+'pg1001/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1001.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1004:
   coadd_x1d,wave,flux,err,path=datdir+'pg1004/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1004.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1116:
   coadd_x1d,wave,flux,err,path=datdir+'pg1116/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1116.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1126:
   coadd_x1d,wave,flux,err,path=datdir+'pg1126/',chan=1,plot=1,bin=3,$
      files=['lcbx01kmq_x1d','lcbx01ktq_x1d','lcbx02hvq_x1d','lcbx02hoq_x1d',$
      'lcbx03bgq_x1d','lcbx03bnq_x1d','lbp408tsq_x1d','lbp408tuq_x1d',$
      'lbp408twq_x1d','lbp408tyq_x1d','lcn701kmq_x1d', 'lcn701koq_x1d']
   openw,lun,specdir+'pg1126.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   coadd_x1d,wave,flux,err,path=datdir+'pg1126/',chan=1,plot=1,bin=3,$
      files=['lcbx01l0q_x1d','lcbx01laq_x1d','lcbx01lhq_x1d','lcbx01ljq_x1d',$
      'lcbx02hxq_x1d','lcbx02i0q_x1d','lcbx02i7q_x1d','lcbx02i9q_x1d',$
      'lcbx03bpq_x1d','lcbx03buq_x1d','lcbx03bxq_x1d','lcbx03c1q_x1d',$
      'lcn701ksq_x1d','lcn701kwq_x1d','lcn701kyq_x1d','lcn701l0q_x1d']
   openw,lun,specdir+'pg1126_cenwave1055.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun

   ;PG1202: Looking at Kriss 2019 FUVA-only data. One of the exposures
   ;(le2j04bwq) has no signal, and the other data are lower-S/N and honestly look
   ;different in the lines from the 2014 data. So given that they are
   ;difficult to incorporate anyway, just going to ignore. If I change my
   ;mind in the future, can do FUVA-only, then FUVB-only, and then just
   ;stitch the text files together.
   ;coadd_x1d,wave,flux,err,path=datdir+'pg1202/',chan=1,plot=1,bin=3,/aonly
   ;openw,lun,specdir+'pg1202_fuva.txt',/get_lun
   ;for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   ;free_lun,lun
   
   ;Then move the FUVA-only files back to their own subdir and run:
   ;coadd_x1d,wave,flux,err,path=datdir+'pg1202/',chan=1,plot=1,bin=3
   ;openw,lun,specdir+'pg1202.txt',/get_lun
   ;for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   ;free_lun,lun
   
   ;PG1211:
   coadd_x1d,wave,flux,err,path=datdir+'pg1211/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1211.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun

   ;PG1226:
   coadd_x1d,wave,flux,err,path=datdir+'pg1226/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1226.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1229:
   coadd_x1d,wave,flux,err,path=datdir+'pg1229/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1229.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;MRK231:
   coadd_x1d,wave,flux,err,path=datdir+'mrk231/',chan=1,plot=1,bin=3
   openw,lun,specdir+'mrk231.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1302:
   coadd_x1d,wave,flux,err,path=datdir+'pg1302/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1302.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1307:
   coadd_x1d,wave,flux,err,path=datdir+'pg1307/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1307.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1309:
   coadd_x1d,wave,flux,err,path=datdir+'pg1309/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1309.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1351:
   coadd_x1d,wave,flux,err,path=datdir+'pg1351/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1351.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1411:
   ;G130M, CENWAVE=1300
   coadd_x1d,wave,flux,err,path=datdir+'pg1411/',chan=1,plot=1,bin=3,$
      files=['lbp414niq_x1d','lbp414nkq_x1d','lbp414nmq_x1d','lbp414noq_x1d',$
      'lc9904a4q_x1d','lc9904abq_x1d','lc9904adq_x1d','lc9904afq_x1d',$
      'ld2n01nlq_x1d','ld2n01nsq_x1d','ld2n01nuq_x1d','ld2n01nwq_x1d',$
      'lddy01hbq_x1d','lddy01hdq_x1d','lddy01hgq_x1d','lddy01hiq_x1d']
   openw,lun,specdir+'pg1411.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   ;G130M, CENWAVE=1096
   ;Added scaling and wavelength xcor:
   coadd_x1d,wave,flux,err,path=datdir+'pg1411/',chan=1,plot=1,bin=3,$
      files=['lc9904akq_x1d','lc9904asq_x1d','lc9904avq_x1d','lc9904b2q_x1d',$
      'ld2n01nzq_x1d','ld2n01o5q_x1d','ld2n01o7q_x1d','ld2n01o9q_x1d',$
      'lddy01hkq_x1d','lddy01hmq_x1d','lddy01hoq_x1d','lddy01hqq_x1d']
   openw,lun,specdir+'pg1411_cenwave1096.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1435:
   coadd_x1d,wave,flux,err,path=datdir+'pg1435/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1435.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1440:
   coadd_x1d,wave,flux,err,path=datdir+'pg1440/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1440.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1448:
   coadd_x1d,wave,flux,err,path=datdir+'pg1448/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1448.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1501:
   coadd_x1d,wave,flux,err,path=datdir+'pg1501/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1501.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1613:
   coadd_x1d,wave,flux,err,path=datdir+'pg1613/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1613.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1617:
   coadd_x1d,wave,flux,err,path=datdir+'pg1617/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1617.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG1626:
   coadd_x1d,wave,flux,err,path=datdir+'pg1626/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg1626.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG2130:
   coadd_x1d,wave,flux,err,path=datdir+'pg2130/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg2130.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG2214:
   coadd_x1d,wave,flux,err,path=datdir+'pg2214/both/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg2214.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG2233:
   coadd_x1d,wave,flux,err,path=datdir+'pg2233/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg2233.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun
   
   ;PG2349:
   coadd_x1d,wave,flux,err,path=datdir+'pg2349/',chan=1,plot=1,bin=3
   openw,lun,specdir+'pg2349.txt',/get_lun
   for i=0,n_elements(wave)-1 do printf,lun,wave[i],flux[i],err[i],format=fmt
   free_lun,lun

end
