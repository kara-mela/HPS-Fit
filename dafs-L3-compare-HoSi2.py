"""
data manipulation (adding, norming, ...) and plotting
substance: Ho2PdSi3
edge: L3 (8071eV, X-Ray Data Booklet)
XAFS sim_col: 1
"""

import pylab as pl
import evaluationtools as et
from kara_tools import TUBAF
from kara_tools import functions as kf
from matplotlib.ticker import FixedLocator, MultipleLocator
from evaluationtools import absorption as ab
# import pyFDMNES as pF

def R_factor(exp, sim, weights=1):
    """
    reliability factor
    Zanazzi1977
    http://journals.iucr.org/q/issues/1965/03/00/a04554/a04554.pdf
    
    How to apply weights inline?
    R**2 = (sum w_i(|F_i^o|-|F_i^c|)**2) / (sum w_i*|F_i^o|**2)
    """
    a = pl.multiply(weights, (abs(exp)-abs(sim))**2)
    b = pl.multiply(weights, abs(exp)**2)
    R = pl.sqrt(sum(a) / sum(b))
    return R

ps = 'TUBA' # plotstyle
pl.matplotlib.rc('font', **{'size':14})

density = 7.6062541648443736
edge = 8071 + 6.3
cut = 42
E_lim = slice(0, 350) # only L3 edge
fact = 1.5 # shift of graphs on y axis

Reflections = {"sat" : "-215", 
               "110" : "220", 
               "001" : "008",
               "301" : "608"}

ExpFunc = {} # Experimental Data as Functions
Sim = {} # Simulated Data

k = {"sat" : 0*fact, 
     "110" : 1*fact, 
     "001" : 2*fact,
     "301" : 3*fact}

# load data
print("loading data...")
for R in Reflections:
    ExpFunc[R] = kf.get_exp(R, norm="mean", crop=E_lim)
    Sim.update(kf.get_sim(R, Reflections, edge))

dE = kf.get_dE(Sim.keys())

# Daten beschneiden und zusammenfuehren
for key in Sim.keys():
    R = key.split('_')[-1]
    if not "FDM" in key:
        continue
    Sim[key][0] += -dE[key] # Energy correction
    if "mod" in key:
        continue
    
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

theta = { # at E = 8075eV
         "sat" : 12.9625957428644,
         "110" : 22.2779615832629, 
         "001" : 11.0648255327577,
         "301" : 43.1643687250282}

# Fitten
g0 = {"sat" : 0.5, 
      "110" : 0.5, 
      "001" : 0.5,
      "301" : 0.01}
fit_para, fit, fitE = {}, {}, {}
for key in Sim:
    R = key.split("_")[-1]
    E, Isim, Abs = Sim[key]
    
    p0 = dict(m=0., n=1., theta=theta[R]/180.*pl.pi, phi=theta[R]/180.*pl.pi, d=pl.inf, 
          Exp=ExpFunc[R], mu=Abs, dE=dE[key], g=g0[R], Isim=Isim)
    myvars = ["m", "n", "theta", "phi", "g"]
    # myvars = ["m", "n", "theta", "phi"]
    fit_para[key] = et.fitls(E, pl.zeros(len(E)), kf.Ext_fit, p0, 
                              myvars, fitalg='simplex', maxfun=1e6, maxiter=1e6)
                              
    fit[key] = kf.Ext_fit(E, diff=False, **fit_para[key].popt)
    fitE[key] = E

kf.make_fit_dat(fit_para, name='dafs', edge='L')



R_fact={}
for key in fit.keys():
    R = key.split('_')[-1]
    # weights = et.gaussian(x=fitE[key], x0=edge, amp=edge, w=1, y0=0.)
    w = 20
    weights = (fitE[key]>(edge-w)) * (fitE[key]<(edge+w))
    R_fact[key] = R_factor(exp=ExpFunc[R](fitE[key]), sim=fit[key], weights=weights)
    print key, R_fact[key]
    
    #ordentliche Ausgabe!!!




        
#----------------------------------------------------------
# Plot fit results
print("plotting...")
f = pl.figure(figsize=(6,10))
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
    
    # plot simulations
    if ("sat" in key and "A" in key) or ("sat" in key and "HS" in key):
        pass
    elif "c1" in key or "c2" in key:
        pl.plot(fitE[key], fit[key]+k[R], 
                  lw=TUBAF.width(ps), color=color)[0]
    else:
        lines[label] = pl.plot(fitE[key], fit[key]+k[R], 
                  lw=TUBAF.width(ps), color=color)[0]
    
for R in ExpFunc:
    # plot experiment
    lines["Experiment"] = pl.plot(ExpFunc[R].x, 
                                  ExpFunc[R].y+k[R], marker='.', color='black')[0]

pl.ylim([-0.1,6.7])
pl.xlim([8040,8160])

my_legend = pl.legend(lines.values(), lines.keys(), 
                      bbox_to_anchor=(1., .81), ncol=1, prop={'size':14})
# pl.legend(lines.values(), lines.keys(), loc=1, prop={'size':14})

# distances of ticks on axes
pl.axes().xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))
pl.axes().yaxis.set_major_locator(MultipleLocator(1.5))

# labels
pl.xlabel('Energy (eV)', fontsize=16)
pl.ylabel('Intensity (a. u.)', fontsize=16)
for R in Reflections:
    pl.text(8150, 0.95+k[R], R)

# border line FDM--Green
pl.plot([edge+cut,edge+cut], [-1, 105], color='.75', lw=TUBAF.width(ps), linestyle='-.')
pl.text(edge+cut+1.5, .07, 'Green', fontsize=14, color='.75')
pl.arrow(edge+cut+2, .02, 8, 0., head_width=0.04, head_length=4, fc='.75', ec='.75')
pl.text(edge+cut-11, .07, 'FDM', fontsize=14, color='.75')
pl.arrow(edge+cut-2, .02, -6, 0., head_width=0.04, head_length=4, fc='.75', ec='.75')
 
pl.savefig('dafs-compare-HoSi2-' + TUBAF.name(ps) + '.pdf', transparent=True)
pl.show()
