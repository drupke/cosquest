; 2016jul14  DSNR  created
pro cos_runplots

   table = '/Users/drupke/Box Sync/cosquest/cosquest.csv'
   plotdir = '/Users/drupke/Box Sync/cosquest/plots/paper/'

   cos_uvlineplots,table,plotdir,'pg0844',/NV
   cos_uvlineplots,table,plotdir,'pg0923',/OVI
   cos_uvlineplots,table,plotdir,'pg0953',/OVI
   cos_uvlineplots,table,plotdir,'pg1004',/OVI
   cos_uvlineplots,table,plotdir,'pg1116',/OVI
   cos_uvlineplots,table,plotdir,'pg1126',/NV
   cos_uvlineplots,table,plotdir,'pg1307',/OVI
   cos_uvlineplots,table,plotdir,'pg1351',/NV
   cos_uvlineplots,table,plotdir,'pg1411',/NV
   cos_uvlineplots,table,plotdir,'pg1440',/NV
   cos_uvlineplots,table,plotdir,'pg1613',/NV,/OVI
   cos_uvlineplots,table,plotdir,'pg1617',/NV,/OVI
   cos_uvlineplots,table,plotdir,'pg2130',/NV
   cos_uvlineplots,table,plotdir,'pg2214',/NV

end
