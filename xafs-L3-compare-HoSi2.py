"""
data manipulation (adding, norming, ...) and plotting
substance: Ho2PdSi3
edge: L3 (8071eV, X-Ray Data Booklet)
XAFS sim_col: 1
"""

from scipy.optimize import curve_fit, leastsq, fmin
from scipy.interpolate import interp1d
import pylab as pl
from itertools import product
from kara_tools import TUBAF
import evaluationtools as et

ps = 'TUBAF'
simplex = False

from matplotlib import rc
rc('font', **{'size':16})

def Icorr(E, m, n, I):
    return (m*(E-E[0]) + n) * I

def residuals(v, Sim, Exp, Energy):
    m, n, dE = v
    if simplex:
        return (((Icorr(Energy, m, n, Sim) - Exp(Energy - dE))**2).sum())
    else:
        return Icorr(Energy, m, n, Sim) - Exp(Energy - dE)

def sim_cut(key, edge, cut, dE):
    # return the index of the biggest energy that is still smaller, than edge + cut
    return max(( idx for idx in range(len(energy[key])) if (energy[key][idx] < edge + cut + dE)), key=lambda idx: idx)

edge = 8071
cut = 45

data, energy, xafs, fit_para, fit = {}, {}, {}, {}, {}

# load data
data["mod"] = pl.loadtxt("modulated-L23-conv_out_conv.txt", skiprows=1)
data["D1_Green"] = pl.loadtxt("D1-L23-Green-conv_out_conv.txt", skiprows=1)
data["D1_FDM"] = pl.loadtxt("MN-v11_conv.txt", skiprows=1)
data["A_Green"] = pl.loadtxt("A-L23-Green-conv_out_conv.txt", skiprows=1)
data["A_FDM"] = pl.loadtxt("A-L23-new-all-conv_out_conv.txt", skiprows=1)
data["HS_Green"] = pl.loadtxt("HoSi2-Green-conv_out_conv.txt", skiprows=1)
data["HS_FDM"] = pl.loadtxt("HoSi2-conv_out_conv.txt", skiprows=1)
exp_data = pl.loadtxt("dafs_hps_sat_psi1_t1_corr.dat", skiprows=1)
# exp_data = pl.loadtxt("dafs_hps_110_ho_r1_corr_enec.dat", skiprows=1)

# energy
for key in data:
    energy[key] = data[key][:,0] + edge
    xafs[key] = data[key][:,1]
en_exp = exp_data[:,0]
Exp = interp1d(en_exp, exp_data[:,-1], kind='linear')

# fit
p0 = [0., 1., 0.2]
for key in xafs:
    args = (xafs[key], Exp, energy[key])
    if simplex:
        fit_para[key] = fmin(residuals, p0, args=args)
    else:
        fit_para[key] = leastsq(residuals, p0, args=args, full_output=True)[0]
    
print "m, n, dE"
header, content = [], []
for key in fit_para:
    print fit_para[key]
    m, n, dE = fit_para[key]
    fit[key] = Icorr(energy[key], m, n, xafs[key])
    header.append(key)
    content.append(dE)
data = dict(zip(header, content))
et.savedat('fit-para.dat', data, xcol=header[0])

# norming
for key in fit: 
    # fit[key] /= max(fit[key])
    fit[key] = (fit[key] - min(fit[key]))/(max(fit[key]) - min(fit[key]))
exp_norm = (Exp(en_exp) - Exp(en_exp[:284]).min())/(Exp(en_exp[:284]).max() - Exp(en_exp[:284]).min())

# cut
idx = {}
models = ['D1', 'A', 'HS']
for i in range(len(models)):
    key_G = models[i] + '_Green'
    key_F = models[i] + '_FDM'
    
    idx_G = sim_cut(key_G, edge, cut, fit_para[key_G][2])
    idx_F = sim_cut(key_F, edge, cut, fit_para[key_F][2])
    # binary arrays indicating wether or not the energy is bigger than cut-energy
    idx[key_G] = (pl.array(range(len(energy[key_G]))) >= idx_G)
    idx[key_F] = -(pl.array(range(len(energy[key_F]))) >= idx_F)
    
    ratio = fit[key_F][idx_G] / fit[key_G][idx_G] # ratio of intensities at cut-energy
    fit[key_G] = fit[key_G]*ratio*idx[key_G] + fit[key_F]*idx[key_F] # combining of both models to one curve
    
    # repairing zeros from last step at cut-energy
    for j in range(len(fit[key_G])):
        if fit[key_G][j] == 0.:
            fit[key_G][j] = fit[key_G][j+1]
idx["mod"] = 1.

# Plot fit results
f = pl.figure()
for key in fit: 
    if "FDM" not in key:
        print key, ps
        if "mod" in key:
            color = TUBAF.gruen(ps)
        elif "D1" in key:
            color = TUBAF.rot(ps)
        elif "A" in key:
            color = TUBAF.orange(ps)
        elif "HS" in key:
            color = TUBAF.blau(ps)
    
        label = key.split('_')[0]

        pl.plot(energy[key] - fit_para[key][2], fit[key], label=label, lw=2*TUBAF.width(ps), color=color)
pl.plot(en_exp, exp_norm, label='Experiment', color='black', marker='.')

pl.ylim([-0.05,1.07])
pl.xlim([8025,8199])
pl.legend(loc=2, prop={'size':12})
pl.xlabel('Energy [eV]', fontsize=18)
pl.ylabel('Intensity [a. u.]', fontsize=18)

# border line FDM--Green
pl.plot([edge+cut,edge+cut], [-1, 105], color='gray', lw=2*TUBAF.width(ps), linestyle='--')
pl.text(edge+cut+5.5, 1., 'Green', fontsize=16, color='0.33')
pl.arrow(edge+cut+6, .97, 15, 0., head_width=0.02, head_length=5, fc='0.33', ec='0.33')
pl.text(edge+cut-20, 1., 'FDM', fontsize=16, color='0.33')
pl.arrow(edge+cut-6, .97, -15, 0., head_width=0.02, head_length=5, fc='0.33', ec='0.33')
 
pl.savefig('xafs-compare-HoSi2-' + TUBAF.name(ps) + '.pdf', transparent=True)
pl.show()

