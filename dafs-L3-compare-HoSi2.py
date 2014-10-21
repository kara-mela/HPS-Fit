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
from kara_tools import TUBAF
import matplotlib
import os
import evaluationtools as et

ps = 'TUBAF'
simplex = True# False#

from matplotlib import rc
rc('font', **{'size':14})

def Icorr(E, Isim, Exp, m=0, n=1, c=0., Abs=1, dE=0, diff=True):
    return (m*E + n) * Isim / Abs**c  - Exp(E - dE) * diff

def sim_cut(key, edge, cut, dE, shift):
    # return the index of the biggest energy that is still smaller, than edge + cut
    return max(( idx for idx in range(len(energy[key])) if (energy[key][idx] < edge + cut + dE - shift)), key=lambda idx: idx)

edge = 8071
cut = 45

reflections = ("301", "sat", "110", "001")
models      = ("mod", "D1", "A", "HS")
algs        = ("FDM", "Green")
  
# Cols = collections.namedtuple("Col", "ss sp")
# cols = dict({"301" : Cols(26,32),
             # # "sat" : Cols(16,23), # only for D1_FDM_sat
             # "sat" : Cols(38,44), # I statt Ic, Abs = 44
             # "110" : Cols(5,11),
             # "001" : Cols(17,23)})

fit_para, fit, exp_norm = {}, {}, {}

# load data
print("loading data...")
all_files = os.listdir(os.curdir)

exp_data = {}
for key in reflections:
    fname = filter(lambda file: key in file and "corr" in file, all_files)[0]
    dummy = et.loaddat(fname)
    exp_data[key] = dict(zip(dummy[1].split(), dummy[0]))


name_dict = dict({"mod_Green"  : "modulated-L23-conv_out_conv.txt",
                  "D1_Green"   : "D1-L23-Green-conv_out_conv.txt",
                  "D1_FDM_sat" : "MN-v16_conv.txt",
                  "D1_FDM"     : "MN-v11_conv.txt",
                  "A_Green"    : "A-L23-Green-conv_out_conv.txt",
                  "A_FDM"      : "A-L23-new-all-conv_out_conv.txt",
                  "HS_Green"   : "HoSi2-Green-conv_out_conv.txt",
                  "HS_FDM"     : "HoSi2-conv_out_conv.txt"})
data = {}
for key in name_dict:
    dummy = et.loaddat(name_dict[key])
    data[key] = dict(zip(dummy[1].split(), dummy[0]))

myvars = ["m", "n", "c", "dE"]
 
# energy/intensity; pre-norming by mean
energy, dafs, Abs = {}, {}, {}
for key, Ref in product(data, reflections):
    if len(key.split('_')) == 2:
        new_key = key + '_' + Ref
    else:
        new_key = key
        Ref = key.split('_')[2]
    
    energy[new_key] = data[key]["Energy"] + edge
    
    I_ss  = "I(" + Ref + ")ss_0"
    I_sp  = "I(" + Ref + ")sp_0"
    Ic_ss = "Ic(" + Ref + ")ss_0"
    Ic_sp = "Ic(" + Ref + ")sp_0"
    A     = "A(" + Ref + ")in_0"
    
    if Ref == "sat":
        dafs[new_key] = data[key][I_ss] + data[key][I_sp]
        Abs[new_key]  = data[key][A]
    else:
        dafs[new_key] = data[key][Ic_ss] + data[key][Ic_sp]
        Abs[new_key] = 1.
    dafs[new_key] /= dafs[new_key].mean()
 
""" 
en_exp, Exp = {}, {}
for key in exp:
    # en_exp[key] = exp[key][:,0]
    # dummy_dafs = exp[key][:,-3] / exp[key][:,-3].mean()
    # Exp[key] = interp1d(en_exp[key], dummy_dafs, kind='linear')
    
    en_exp[key] = exp_data[key]
  
# fit
Abs = {}
print("fitting...")
for key in dafs:
    if "sat" in key and "D1" in key:
        if "FDM" in key: 
            print 'FDM', key
            Abs[key] = data[key][:,21]
        else: 
            print 'Green', key
            Abs[key] = data[key[:-4]][:,42]
    elif "sat" in key and "mod" in key:
        Abs[key] = 1.
    else:
        Abs[key] = 1.
for key in dafs:
    # p0 = dict(m=0.01, n=1., c=1., dE=0., Exp=Exp[key.split('_')[2]], 
      # Abs=Abs[key], Isim = dafs[key])
    p0 = dict(m=1., n=0., c=1., dE=1.5, Exp=Exp[key.split('_')[2]], 
      Abs=Abs[key], Isim = dafs[key])
    
    name = key.split('_')[-1]
    args = (dafs[key], Exp[name], energy[key], Abs[key])
    
    fit_para[key] = et.fitls(energy[key], pl.zeros(len(energy[key])), Icorr, p0, 
      myvars, fitalg="simplex")
    
    fit[key] = Icorr(energy[key], diff=False, **fit_para[key].popt)

# norming
print("norming...")

# Falsche Normierung fuer DAFS!!!

for key in fit: 
    fit[key] = (fit[key] - min(fit[key]))/(max(fit[key]) - min(fit[key]))
    
    refl = key.split('_')[2]
    f = Exp[refl]
    en = en_exp[refl]
    exp_norm[refl] = (f(en) - f(en[:284]).min())/(f(en[:284]).max() - f(en[:284]).min())
    
k, shift = {}, {}
j = 0
for key in cols:
    k[key] = j
    shift[key] = -23.5 if key == "001" else 0.
    j += 1.

# cut
print("cutting...")
idx = {}
models = ['D1', 'A', 'HS']
for i, j in product(range(len(models)), cols):
    key_G = models[i] + '_Green_' + j
    key_F = models[i] + '_FDM_' + j
    
    # index of cut energy
    idx_G = sim_cut(key_G, edge, cut, fit_para[key_G].popt["dE"], shift[j])
    idx_F = sim_cut(key_F, edge, cut, fit_para[key_F].popt["dE"], shift[j])
    
    # binary arrays indicating wether or not the energy is bigger than cut energy
    idx[key_G] = (pl.array(range(len(energy[key_G]))) >= idx_G)
    idx[key_F] = -(pl.array(range(len(energy[key_F]))) >= idx_F)
    
    # ratio of intensities at cut-energy
    ratio = fit[key_F][idx_G] / fit[key_G][idx_G] 
    
    # create array consisting of FDM data upto cut energy, followed by Green data
    dafs_dummy, energy_dummy = [], []
    for i in range(len(idx[key_F])):
        if idx[key_F][i] != 0:
            dafs_dummy.append(fit[key_F][i])
            energy_dummy.append(energy[key_F][i])
    for i in range(len(idx[key_G])):
        if idx[key_G][i] != 0:
            dafs_dummy.append(fit[key_G][i]*ratio)
            energy_dummy.append(energy[key_G][i])
    fit[key_G] = pl.array(dafs_dummy)
    energy[key_G] = pl.array(energy_dummy)
idx["mod"] = 1.

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
            print key, fit_para[key].popt["dE"], fit_para[key].popt["c"]
            
            pl.plot(energy[key] + fit_para[key].popt["dE"], fit[key] + k[refl], 
              lw=2*TUBAF.width(ps), color=color)
        
    # plot experiment
    pl.plot(en_exp[refl] + shift[refl], exp_norm[refl] + k[refl], color='black', marker='.')

pl.ylim([-0.1,4.3])
pl.xlim([8025,8199])

# virtual curves for creating legend
my_exp = matplotlib.lines.Line2D(en_exp["001"], exp_norm["001"], marker='.', color='black')
my_mod = matplotlib.lines.Line2D(energy["mod_Green_001"], fit["mod_Green_001"], lw=2*TUBAF.width(ps), color=TUBAF.gruen(ps))
my_D1 = matplotlib.lines.Line2D(energy["D1_Green_001"], fit["D1_Green_001"], lw=2*TUBAF.width(ps), color=TUBAF.rot(ps))
my_A = matplotlib.lines.Line2D(energy["A_Green_001"], fit["A_Green_001"], lw=2*TUBAF.width(ps), color=TUBAF.orange(ps))
my_HS = matplotlib.lines.Line2D(energy["HS_Green_001"], fit["HS_Green_001"], lw=2*TUBAF.width(ps), color=TUBAF.blau(ps))

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
for key in cols:
    pl.text(8160, .6+k[key], key)

# border line FDM--Green
pl.plot([edge+cut,edge+cut], [-1, 105], color='gray', lw=2*TUBAF.width(ps), linestyle='--')
pl.text(edge+cut+5.5, 4.15, 'Green', fontsize=16, color='0.33')
pl.arrow(edge+cut+6, 4.12, 25, 0., head_width=0.035, head_length=5, fc='gray', ec='gray')
pl.text(edge+cut-30, 4.15, 'FDM', fontsize=16, color='0.33')
pl.arrow(edge+cut-6, 4.12, -25, 0., head_width=0.035, head_length=5, fc='gray', ec='gray')
 
pl.savefig('dafs-compare-HoSi2-' + TUBAF.name(ps) + '.pdf', transparent=True)
pl.show()
"""
