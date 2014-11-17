"""
data manipulation (adding, norming, ...) and plotting
substance: Ho2PdSi3
edge: L3 (8071eV, X-Ray Data Booklet)
XAFS sim_col: 1
"""

from scipy.optimize import curve_fit, leastsq, fmin
from scipy.interpolate import interp1d
import pylab as pl
import collections
from matplotlib.ticker import MultipleLocator
from itertools import product
from kara_tools import TUBAF, io
import matplotlib
import os
import evaluationtools as et

ps = 'TUBAF'
simplex = False# True#

from matplotlib import rc
rc('font', **{'size':14})

def Icorr(E, Isim, Exp, m=0, n=1, c=0., Abs=1, diff=True):
    return (m*(E-E[0]) + n) * Isim / Abs**c  - Exp(E) * diff

def sim_cut(key, edge, cut, dE, shift):
    # return the index of the biggest energy that is still smaller, than edge + cut
    return max(( idx for idx in range(len(Energy[key])) if (Energy[key][idx] < edge + cut + dE[key] - shift)), key=lambda idx: idx)

def get_dE(keys):
    dE_xafs = et.loaddat('fit-para.dat', todict=True, comment='#')
    dE = {}
    for idx, key in product(dE_xafs, keys):
        if idx.split('_')[0] in key and idx.split('_')[-1] in key:
            dE[key] = dE_xafs[idx][0]
    return dE

def get_exp(R, DIR = os.curdir):
    flist = os.listdir(DIR)
    fname = filter(lambda s: s.startswith("dafs_hps_%s"%R), flist)
    assert len(fname) == 1, "More than 1 file found."
    fname = fname[0]
    data = et.loaddat(fname, todict=True, comment='')
    return interp1d(data["Energy"], data["bragg"], kind='linear')
    
def idx_cut_energy(edge, cut, energy, shift, dE):
    """
    determing closest index to cut energy, concerning shift and dE
    """
    return list(pl.round_(energy)).index(edge + cut - shift + dE)

def patching(fit_F, fit_G, idx_F, idx_G, e_F, e_G):
    """
    ratio for FDM AND Green sim -> adapting to exp
    """
    ratio = fit_F[idx_F] / fit_G[idx_G] 
    
    dafs_dummy, energy_dummy = [], []
    
    for i in range(idx_F):
        dafs_dummy.append(fit_F[i])
        energy_dummy.append(e_F[i])
    
    for i in range(len(e_G) - idx_G):
        dafs_dummy.append(fit_G[i]*ratio)
        energy_dummy.append(e_G[i])
        
    return pl.array(dafs_dummy), pl.array(energy_dummy)

    
edge = 8071
cut = 45

# myvars = ["m", "n", "c", "dE"]
myvars = ["m", "n", "c"]
  
fit_para, fit, exp_norm = {}, {}, {}

# load data
Reflections = {"301" : "608", 
               "sat" : "-215", 
               "110" : "220", 
               "001" : "008"}

Models   = {'A-L23-Green-conv_out_conv.txt'  : 'A_Green', 
            'A-L23-new-all-conv_out_conv.txt': 'A_FDM', 
            'D1-L23-Green-conv_out_conv.txt' : 'D1_Green', 
            'HoSi2-Green-conv_out_conv.txt'  : 'HS_Green', 
            'HoSi2-conv_out_conv.txt'        : 'HS_FDM', 
            'MN-v11_conv.txt'                : 'D1_FDM', 
            'MN-v16_conv.txt'                : 'D1_FDM', 
            'modulated-L23-conv_out_conv.txt': 'mod_Green'}

ExpFunc = {}
Energy = {}
Sim = {}
Abs = {}

def get_sim(R, Reflections, simfile):
    useabs = R=="sat"
    Ref = Reflections[R] if not ("HS_" in Models[simfile]) else R
    
    try:
        data = io.FDMNES.loadDAFS(simfile, Ref, absorption=useabs)
        if useabs:
            Abs[key] = data[2]
        return data[1] / data[1].mean(), data[0] + edge
    except ValueError: #not found in this model?
        print("Reflection %s not found in file %s"%(R, simfile))

print("loading data...")
for R in Reflections:
    ExpFunc[R] = get_exp(R)

    for simfile in Models:
        key = "_".join([Models[simfile], R])
        
        # print key
        # Sim[key], Energy[key] = get_sim(R, Reflections, simfile)
        
        useabs = R=="sat"
        Ref = Reflections[R] if not ("HS_" in Models[simfile]) else R
        
        try:
            data = io.FDMNES.loadDAFS(simfile, Ref, absorption=useabs)
        except ValueError: #not found in this model?
            print("Reflection %s not found in file %s"%(R, simfile))
            continue
        
        if useabs:
            Abs[key] = data[2]
        
        Sim[key] = data[1] / data[1].mean()
        Energy[key] = data[0] + edge

print "\ndone\n"

dE = get_dE(Sim.keys())  
  
# fitting
for key in Sim:
    R = key.split("_")[-1]
    
    p0 = dict(m=0.01, n=1.01, c=1.01, Exp=ExpFunc[R], Isim = Sim[key])
    
    if key in Abs:
        p0["Abs"] = Abs[key]
    
    fit_para[key] = et.fitls(Energy[key] - dE[key], pl.zeros(len(Energy[key])), Icorr,
      p0, myvars, fitalg="simplex")
    
    # print key, fit_para[key].popt["c"], fit_para[key].popt["m"], fit_para[key].popt["n"]
    
    fit[key] = Icorr(Energy[key], diff=False, **fit_para[key].popt)

# norming
print("norming...")
for key in fit: 
    fit[key] /= fit[key].max()
    
    R = key.split('_')[-1]
    f = ExpFunc[R]
    # exp_norm[R] = f(ExpEn[R])/f(ExpEn[R][:284]).max()
    exp_norm[R] = interp1d(f.x[:284], f.y[:284], kind='linear')
    
k, shift = {}, {}
j = 0
for key in Reflections:
    k[key] = j
    shift[key] = -23.5 if key == "001" else 0.
    j += 1.

# cut
print("cutting...")
for key_F in Sim:
    if not "FDM" in key:
        continue
    key_G = key_F.replace("FDM", "Green")
    
    idx_F = idx_cut_energy(edge, cut, Energy[key_F], shift[key_F], dE[key_F])
    idx_G = idx_cut_energy(edge, cut, Energy[key_G], shift[key_G], dE[key_G])
    
    fit[key_G], Energy[key_G] = patching(fit[key_F], fit[key_G], 
      idx_F, idx_G, Energy[key_F], Energy[key_G])

#----------------------------------------------------------
# Plot fit results
print("plotting...")
f = pl.figure(figsize=(5,10))
for key in fit:
    refl = key.split('_')[2]
    if "FDM" not in key:
        if "mod" in key:
            color = TUBAF.gruen(ps)
        elif "D1" in key:
            color = TUBAF.rot(ps)
        elif "A" in key:
            color = TUBAF.orange(ps)
        elif "HS" in key:
            color = TUBAF.blau(ps)
        
        # plot simulations
        if ("sat" in key and "A" in key) or ("sat" in key and "HS" in key):
            None
        else:
            pl.plot(Energy[key]-dE[key], fit[key] + k[refl], lw=2*TUBAF.width(ps), color=color)
        
    # plot experiment
    print refl
    pl.plot(Energy[key] + shift[refl], exp_norm[refl](Energy[key]) + k[refl], color='black', marker='.')

# test
pl.ylim([-0.1,4.3])
pl.xlim([8025,8199])

# virtual curves for creating legend
my_exp = matplotlib.lines.Line2D(*Exp["001"], marker='.', color='black')
my_mod = matplotlib.lines.Line2D(Energy["mod_Green_001"], fit["mod_Green_001"], lw=2*TUBAF.width(ps), color=TUBAF.gruen(ps))
my_D1 = matplotlib.lines.Line2D(Energy["D1_Green_001"], fit["D1_Green_001"], lw=2*TUBAF.width(ps), color=TUBAF.rot(ps))
my_A = matplotlib.lines.Line2D(Energy["A_Green_001"], fit["A_Green_001"], lw=2*TUBAF.width(ps), color=TUBAF.orange(ps))
my_HS = matplotlib.lines.Line2D(Energy["HS_Green_001"], fit["HS_Green_001"], lw=2*TUBAF.width(ps), color=TUBAF.blau(ps))

my_legend = pl.legend((my_exp, my_mod, my_D1, my_A, my_HS), 
  (r"Experiment", r"model mod", r"model $D1$", r"model $A$", r"HoSi$_2$"), 
  bbox_to_anchor=(0., 1.005, 1., .065), loc=3,
  ncol=2, mode="expand", borderaxespad=0., prop={'size':14})

# distances of ticks on axes
pl.axes().xaxis.set_major_locator(MultipleLocator(50))
pl.axes().xaxis.set_minor_locator(MultipleLocator(25))
pl.axes().yaxis.set_minor_locator(MultipleLocator(0.5))

# labels
pl.xlabel('Energy [eV]', fontsize=18)
pl.ylabel('Intensity [a. u.]', fontsize=18)
for key in Reflections:
    pl.text(8160, .6+k[key], key)

# border line FDM--Green
pl.plot([edge+cut,edge+cut], [-1, 105], color='gray', lw=2*TUBAF.width(ps), linestyle='--')
pl.text(edge+cut+5.5, 4.15, 'Green', fontsize=16, color='0.33')
pl.arrow(edge+cut+6, 4.12, 25, 0., head_width=0.035, head_length=5, fc='gray', ec='gray')
pl.text(edge+cut-30, 4.15, 'FDM', fontsize=16, color='0.33')
pl.arrow(edge+cut-6, 4.12, -25, 0., head_width=0.035, head_length=5, fc='gray', ec='gray')
 
pl.savefig('dafs-compare-HoSi2-' + TUBAF.name(ps) + '.pdf', transparent=True)
pl.show()
