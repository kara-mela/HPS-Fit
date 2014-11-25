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
from functions import *

MultipleLocator = pl.matplotlib.ticker.MultipleLocator

ps = 'TUBAF'
pl.matplotlib.rc('font', **{'size':14})

def Icorr(E, Isim, Exp, dE=0, m=0, n=1, c=0., Abs=1, diff=True):
    return (m*(E-E[0]) + n) * Isim / Abs**c  - Exp(E - dE) * diff
    
edge = 8071
cut = 45
E_lim = slice(0, 350) # only L3 edge

#myvars = ["m", "n", "c"]
myvars = ["n", "m"]

shift = collections.defaultdict(int)
shift["001"] = -23.5

Reflections = {"301" : "608", 
               "sat" : "-215", 
               "110" : "220", 
               "001" : "008"}

ExpFunc = {} # Experimental Data as Functions
Sim = {} # Simulated Data

k = dict(zip(Reflections.keys(), 2*pl.arange(len(Reflections))))


# load data
print("loading data...")
for R in Reflections:
    ExpFunc[R] = get_exp(R, norm="mean", crop=E_lim)
    Sim.update(get_sim(R, Reflections, edge))

dE = get_dE(Sim.keys())

for key in Sim:
    Sim[key][0] += dE[key] # Energy correction
# Daten beschneiden und zusammenfuehren
idx = {}
for key in Sim.keys():
    if not "FDM" in key:
        continue
    Sim[key][0] += dE[key] # Energy correction
    if "mod" in key:
        continue
        #limits = slice(0,-1)
    R = key.split('_')[-1]
    
    keyG = key.replace("FDM", "Green")
    idxF = Sim[key][0] <= (edge + cut)
    idxG = Sim[keyG][0] > (edge + cut)
    for i in [1,2]:
        ratio = Sim[keyG][i,idxG][0] / Sim[key][i,~idxF][0]
        Sim[keyG][i,idxG] /= ratio
    Sim[key] = pl.hstack((Sim[key][:,idxF], Sim[keyG][:,idxG]))
    Sim.pop(keyG)

for key in Sim.keys():
    if Sim[key][1].max() < 1e-10:
        Sim.pop(key) # No Intensity
    else:
        Sim[key][1] /= Sim[key][1].mean() # normalize
        Sim[key][2] /= Sim[key][2].mean()


# Fitten
fit_para, fit = {}, {}
fitE = {}
for key in Sim:
    R = key.split("_")[-1]
    E, Isim, Abs = Sim[key]
    p0 = dict(m=0.0, n=1., c=1., Exp=ExpFunc[R], Isim=Isim, Abs=Abs)
    
    fit_para[key] = et.fitls(E, pl.zeros(len(E)), Icorr, p0, myvars, 
                             fitalg="simplex")
    
    # if "sat" in key and "D1" in key:
        # print "   sat-Parameter:", fit_para[key].popt["m"], 
          # fit_para[key].popt["n"], fit_para[key].popt["c"]
    
    fit[key] = Icorr(E, diff=False, **fit_para[key].popt)
    fitE[key] = E

#----------------------------------------------------------
# Plot fit results
print("plotting...")
f = pl.figure(figsize=(5,10))
lines = {}
for key in fit:
    R = key.split('_')[2]
    if "mod" in key:
        color = TUBAF.gruen(ps)
        label = "model mod"
    elif "D1" in key:
        color = TUBAF.rot(ps)
        label = "model D1"
    elif "A" in key:
        color = TUBAF.orange(ps)
        label = "model A"
    elif "HS" in key:
        color = TUBAF.blau(ps)
        label = "HoSi$_2$"
    
    # plot simulations
    if ("sat" in key and "A" in key) or ("sat" in key and "HS" in key):
        pass
    else:
        lines[label] = pl.plot(fitE[key], fit[key]+k[R], 
                               lw=2*TUBAF.width(ps), color=color)[0]
    
for R in ExpFunc:
    # plot experiment
    lines["Experiment"] = pl.plot(ExpFunc[R].x+shift[R], 
                                  ExpFunc[R].y+k[R], '.k')[0]

# test
pl.ylim([-0.1,2*4.3])
pl.xlim([8025,8199])

my_legend = pl.legend(lines.values(), lines.keys(), 
                      bbox_to_anchor=(0., 1.005, 1., .065), loc=3, ncol=2, 
                      mode="expand", borderaxespad=0., prop={'size':14})

# distances of ticks on axes
pl.axes().xaxis.set_major_locator(MultipleLocator(50))
pl.axes().xaxis.set_minor_locator(MultipleLocator(25))
pl.axes().yaxis.set_minor_locator(MultipleLocator(0.5))

# labels
pl.xlabel('Energy [eV]', fontsize=18)
pl.ylabel('Intensity [a. u.]', fontsize=18)
for R in Reflections:
    pl.text(8175, 0.95+k[R], R)

# border line FDM--Green
pl.plot([edge+cut,edge+cut], [-1, 105], color='gray', lw=2*TUBAF.width(ps), linestyle='--')
pl.text(edge+cut+5.5, 2*4.15, 'Green', fontsize=16, color='0.33')
pl.arrow(edge+cut+6, 2*4.12, 25, 0., head_width=0.05, head_length=5, fc='gray', ec='gray')
pl.text(edge+cut-30, 2*4.15, 'FDM', fontsize=16, color='0.33')
pl.arrow(edge+cut-6, 2*4.12, -25, 0., head_width=0.05, head_length=5, fc='gray', ec='gray')
 
pl.savefig('dafs-compare-HoSi2-' + TUBAF.name(ps) + '.pdf', transparent=True)
pl.show()
