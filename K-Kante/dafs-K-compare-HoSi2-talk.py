"""
data manipulation (adding, norming, ...) and plotting
substance: Ho2PdSi3
edge: L3 (8071eV, X-Ray Data Booklet)
XAFS sim_col: 1
"""

import os
import pylab as pl
import collections
import evaluationtools as et
from kara_tools import TUBAF
from kara_tools import functions as f
from matplotlib.ticker import FixedLocator, MultipleLocator

MultipleLocator = pl.matplotlib.ticker.MultipleLocator

ps = 'TUBAF'
pl.matplotlib.rc('font', **{'size':14})

edge = 24365
cut = 42
E_lim = slice(0, -1)
fact=1.

myvars = ["n", "m"]

Reflections = {"001" : "008",
               "301" : "608"}

ExpFunc = {} # Experimental Data as Functions
Sim = {} # Simulated Data

k = dict(zip(Reflections.keys(), fact*pl.arange(len(Reflections))))


# load data
print("loading data...")
for R in Reflections:
    ExpFunc[R] = f.get_exp(R, end="", norm="mean", crop=E_lim)
    Sim.update(f.get_sim(R, Reflections, edge, symbol="K"))

dE = f.get_dE(Sim.keys())

for key in Sim.keys():
    if Sim[key][1].max() < 1e-10:
        print key
        Sim.pop(key) # No Intensity
    else:
        Sim[key][1] /= Sim[key][1].mean() # normalize
        Sim[key][2] /= Sim[key][2].mean()
        
# Fitten
fit_para, fit = {}, {}
fitE = {}
# don't forget the satellite with c=1 and c=2
for key in Sim:
    R = key.split("_")[-1]
    E, Isim, Abs = Sim[key]
    p0 = dict(m=0.0, n=1., c=1., Exp=ExpFunc[R], Isim=Isim, Abs=Abs)

    fit_para[key] = et.fitls(E, pl.zeros(len(E)), f.Icorr, p0, 
                             myvars + ["c"] * (R=="sat"), 
                             # myvars + ["c"] * 1, 
                             fitalg="simplex")
    print fit_para[key].popt["c"]
    
    fit[key] = f.Icorr(E, diff=False, **fit_para[key].popt)
    fitE[key] = E
    
    if ("D1" in key or "mod" in key) and "sat" in key:
        nkey = key + "_c1"
        fit_para[nkey] = fit_para[key]
        fit_para[nkey].popt['c'] = 1.
        fit[nkey] = f.Icorr(E, diff=False, **fit_para[nkey].popt)
        fitE[nkey] = E
        
        nkey = key + "_c2"
        fit_para[nkey] = fit_para[key]
        fit_para[nkey].popt['c'] = 2.
        fit[nkey] = f.Icorr(E, diff=False, **fit_para[nkey].popt)
        fitE[nkey] = E

f.make_fit_dat(fit_para, name='dafs')

#----------------------------------------------------------
# Plot fit results
print("plotting...")
# f = pl.figure(figsize=(7,7))
# f = pl.figure(figsize=(10,20))
lines = {}
for key in fit:
    R = key.split('_')[2]
    if "mod" in key:
        color = TUBAF.color(ps)['g']
        label = "model mod"
    elif "D1" in key:
        color = TUBAF.color(ps)['r']
        label = "model $D_1$"
    elif "A" in key:
        color = TUBAF.color(ps)['o']
        label = "model $A$"
    elif "HS" in key:
        color = TUBAF.color(ps)['b']
        label = "HoSi$_2$"
    if "c1" in key:
        ls = ':'
    elif "c2" in key:
        ls = '--'
    else:
        ls = '-'
    
    # plot simulations
    if ("sat" in key and "A" in key) or ("sat" in key and "HS" in key):
        pass
    elif "c1" in key or "c2" in key:
        pl.plot(fitE[key], fit[key]+k[R], 
                  lw=TUBAF.width(ps), color=color, ls=ls)[0]
    else:
        lines[label] = pl.plot(fitE[key], fit[key]+k[R], 
                  lw=TUBAF.width(ps), color=color, ls=ls)[0]
    
for R in ExpFunc:
    # plot experiment
    lines["Experiment"] = pl.plot(ExpFunc[R].x, 
                                  ExpFunc[R].y+k[R], marker='.', color='black')[0]

# pl.ylim([-0.1,fact*4.3])
pl.ylim([-0.05,2.45])
pl.xlim([24310,24549])
pl.axes().yaxis.set_major_locator(MultipleLocator(1))

pl.legend(lines.values(), lines.keys(), loc=3, prop={'size':12})

# labels
pl.xlabel('Energy (eV)', fontsize=16)
pl.ylabel('Intensity (a. u.)', fontsize=16)
for R in Reflections:
    pl.text(edge + 150, 0.95+k[R], R)
 
pl.savefig('dafs-compare-K-HoSi2-' + TUBAF.name(ps) + '-talk.pdf', transparent=True)
pl.show()
