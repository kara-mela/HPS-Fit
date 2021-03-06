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

MultipleLocator = pl.matplotlib.ticker.MultipleLocator

ps = 'TUBA'
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
    ExpFunc[R] = kf.get_exp(R, end="", norm="mean", crop=E_lim)
    Sim.update(kf.get_sim(R, Reflections, edge, symbol="K"))

dE = kf.get_dE(Sim.keys())

for key in Sim.keys():
    if Sim[key][1].max() < 1e-10:
        print key
        Sim.pop(key) # No Intensity
    else:
        Sim[key][1] /= Sim[key][1].mean() # normalize
        Sim[key][2] /= Sim[key][2].mean()
        
# def Abs_fit(E, Exp, mu, phi, theta, m=0., n=1., d=pl.inf, dE=0., diff=True):
    # """
    # following Booth2005, equation (1)
    # eps_a = 1
    # mu_a = 1
    # mu_t = mu_f, due to elastical scattering
    # """
    # g = pl.sin(phi) / pl.sin(theta)
    # mf = (m*(E-E[0]) + n) * Isim # linear machine function
    # I = mf / (mu + g*mu) * (1 - pl.exp(-d*mu*(1./pl.sin(phi) + 1./pl.sin(theta)))) - Exp(E - dE) * diff
    # return I

theta = { # at E = 24365eV 
         "001" : 3.64679316745243,
         "301" : 13.1040973558876}
    
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
    fit_para[key] = et.fitls(E, pl.zeros(len(E)), kf.Ext_fit, p0, 
                              myvars, fitalg='simplex', maxfun=1e6, maxiter=1e6)
    fit[key] = kf.Ext_fit(E, diff=False, **fit_para[key].popt)
    fitE[key] = E

R_fact = {}
for key in fit.keys():
    R = key.split('_')[-1]
    # weights = et.gaussian(x=fitE[key], x0=edge, amp=edge, w=1, y0=0.)
    w = 20
    weights = (fitE[key]>(edge-w)) * (fitE[key]<(edge+w))
    R_fact[key] = kf.R_factor(exp=ExpFunc[R](fitE[key]), sim=fit[key], weights=weights)

kf.make_fit_dat(fit_para, name='dafs', edge='K', R_fact=R_fact)



#----------------------------------------------------------
# Plot fit results
print("plotting...")
f = pl.figure(figsize=(6,6))
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
# pl.ylim([.55,2.95])
pl.ylim([.55,2.4])
pl.xlim([24310,24549])

my_legend = pl.legend(lines.values(), lines.keys(), 
                      bbox_to_anchor=(1., .57), ncol=1, prop={'size':14})
# pl.legend(lines.values(), lines.keys(), loc=3, prop={'size':14})

# labels
pl.xlabel('Energy (eV)', fontsize=16)
pl.ylabel('Intensity (a. u.)', fontsize=16)
for R in Reflections:
    pl.text(edge + 150, 0.9+k[R], R)
 
pl.savefig('dafs-compare-K-HoSi2-' + TUBAF.name(ps) + '.pdf', transparent=True)
pl.show()
