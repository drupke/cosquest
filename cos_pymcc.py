#!/Users/drupke/opt/anaconda3/envs/astroconda/bin python
from astropy.io import ascii
from astropy.table import Table
import numpy as np
from pymccorrelation import pymccorrelation

Nperturb = 1000

datadir = '/Users/drupke/Box Sync/cosquest/plots/correlations/'

data = ascii.read(datadir+'weq_vs_nhxray_dat.txt')
res = pymccorrelation(data['xdat'],data['ydat'],data['xerr'],data['yerr'],
                      data['xl'],data['yl'],coeff='kendallt',Nperturb=Nperturb,
                      return_dist=True)
# Compute pval from cc distribution
scc = np.sort(res[2])
izero = np.searchsorted(scc,0.)
# case of positive cc, all > 0
# case of negative cc, all < 0
if izero == 0 or izero == scc.size-1:
    pvaldist = 1./float(Nperturb)
    pvalul = True
else:
    pvalul = False
    pvaldist = float(izero)/float(scc.size)
    if res[0][1] < 0.:
        pvaldist = 1. - pvaldist

outtab = Table()
outtab['cc:16/50/84'] = res[0]
outtab['pval:16/50/84'] = res[1]
outtab['pval:calc/ul'] = [pvaldist,pvalul,'N/A']
outtab.write(datadir+'weq_vs_nhxray_stat.txt', format='ascii', delimiter='\t', overwrite=True)


data = ascii.read(datadir+'vwtavg_vs_nhxray_dat.txt')
res = pymccorrelation(data['xdat'],data['ydat'],data['xerr'],data['yerr'],
                      data['xl'],data['yl'],coeff='kendallt',Nperturb=Nperturb,
                      return_dist=True)
# Compute pval from cc distribution
scc = np.sort(res[2])
izero = np.searchsorted(scc,0.)
# case of positive cc, all > 0
# case of negative cc, all < 0
if izero == 0 or izero == scc.size-1:
    pvaldist = 1./float(Nperturb)
    pvalul = True
else:
    pvalul = False
    pvaldist = float(izero)/float(scc.size)
    if res[0][1] < 0.:
        pvaldist = 1. - pvaldist

outtab = Table()
outtab['cc:16/50/84'] = res[0]
outtab['pval:16/50/84'] = res[1]
outtab['pval:calc/ul'] = [pvaldist,pvalul,'N/A']
outtab.write(datadir+'vwtavg_vs_nhxray_stat.txt', format='ascii', delimiter='\t', overwrite=True)



data = ascii.read(datadir+'vwtrms_vs_nhxray_dat.txt')
res = pymccorrelation(data['xdat'],data['ydat'],data['xerr'],data['yerr'],
                      data['xl'],data['yl'],coeff='kendallt',Nperturb=Nperturb,
                      return_dist=True)
# Compute pval from cc distribution
scc = np.sort(res[2])
izero = np.searchsorted(scc,0.)
# case of positive cc, all > 0
# case of negative cc, all < 0
if izero == 0 or izero == scc.size-1:
    pvaldist = 1./float(Nperturb)
    pvalul = True
else:
    pvalul = False
    pvaldist = float(izero)/float(scc.size)
    if res[0][1] < 0.:
        pvaldist = 1. - pvaldist
outtab = Table()
outtab['cc:16/50/84'] = res[0]
outtab['pval:16/50/84'] = res[1]
outtab['pval:calc/ul'] = [pvaldist,pvalul,'N/A']
outtab.write(datadir+'vwtrms_vs_nhxray_stat.txt', format='ascii', delimiter='\t', overwrite=True)
