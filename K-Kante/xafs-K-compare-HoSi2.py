
"""
data manipulation (adding, norming, ...) and plotting
substance: Ho2PdSi3
edge: K (24350eV, X-Ray Data Booklet)
XAFS sim_col: 1
"""

from scipy.optimize import curve_fit, leastsq, fmin
from scipy.interpolate import interp1d
from matplotlib.ticker import FixedLocator, MultipleLocator
import pylab as pl
from itertools import product
import numpy as np
from kara_tools import TUBAF
from kara_tools import functions as f
import evaluationtools as et


def r_value(exp, sim):
    """
    reference: http://reference.iucr.org/dictionary/R_factor
    problem of unequal x-axis: different zeroes, different steps
    """
    r = sum(abs(exp - sim)) / sum(abs(sim))
    return r




ps = 'TUBA'
simplex = False

from matplotlib import rc
rc('font', **{'size':14})

edge = 24365
plot_shift = 16 # shift of features depending on chosen exp

myvars = ["n", "m", "dE"]

data, energy, xafs, fit_para, fit, fitE = {}, {}, {}, {}, {}, {}

# load data
exp_data    = pl.loadtxt('dafs_hps_fluo_pd.dat', skiprows=1)
data['A']   = pl.loadtxt('A-K_conv_out_conv.txt', skiprows=1) 
data['D1']  = pl.loadtxt('D1-K_conv_out_conv.txt', skiprows=1)
# data['mod'] = pl.loadtxt('mod-K_conv_out.txt', skiprows=1)
data['mod'] = pl.loadtxt('mod-K_out_conv.txt', skiprows=1)

# energy
for key in data:
    energy[key] = data[key][:,0] + edge
    xafs[key] = data[key][:,1]
en_exp = exp_data[:,0]
Exp = interp1d(en_exp, exp_data[:,1], kind='linear') # Energy, flutot

for key in xafs.keys():
    if xafs[key].max() < 1e-10:
        xafs.pop(key) # No Intensity
    else:
        xafs[key] /= xafs[key].mean()

# fit
print "m, n, dE"
header, content = [], []
for key in xafs:
    p0 = dict(m=0.0, n=1., c=0., dE=0., Exp=Exp, Isim=xafs[key])
    fit_para[key] = et.fitls(energy[key], pl.zeros(len(energy[key])), 
                         f.Icorr, p0, myvars, fitalg="simplex")
    print fit_para[key].popt["c"]

# for key in fit_para:
    m, n, dE = fit_para[key].popt["m"], fit_para[key].popt["n"], fit_para[key].popt["dE"]
    print m, n, dE
    
    # dE_manu = 2.
    # dE += dE_manu
    # fit_para[key].popt["dE"] += dE_manu
    
    fit[key] = f.Icorr(energy[key], diff=False, **fit_para[key].popt)
    fitE[key] = energy[key]
    
    header.append(key)
    if "mod" in key:
        content.append(dE)
    else:
        content.append(dE)
    
data = dict(zip(header, content))
et.savedat('fit-para.dat', data, xcol=header[0])

# norming
for key in fit: 
    fit[key] = (fit[key] - min(fit[key]))/(np.mean(fit[key][:]) - min(fit[key]))
exp_norm = (Exp(en_exp) - Exp(en_exp[50:]).min())/(np.mean(Exp(en_exp[:450])) - Exp(en_exp[50:]).min())

# Plot fit results
ax1 = pl.axes([.1, .1, .8, .8])

# oscillation labels
def plot_markers(ax):
    # feature markers
    my_labels =             [ '$B$', '$C_1$', '$C_2$', '$C_3$', '$C_4$', '$C_5$']
    my_energies = pl.array( [24.346,  24.360,  24.384,  24.423,  24.480, 24.511])
    my_energies *= 1000
    my_energies += plot_shift
    for i in range(4):
        for line in range(len(my_labels)):  
            if i == 0:
                pl.text(my_energies[line]+.8, 1.825, my_labels[line], fontsize=16)
            
            pl.plot([my_energies[line],my_energies[line]], [-1, 20], 
                     color='gray', lw=0.5*TUBAF.width(ps), linestyle='--')

plot_markers(ax1)

def my_plot(ax, fit):
    for key in fit: 
        if "FDM" not in key:
            if "mod" in key:
                color = TUBAF.color(ps)['g']
                label = 'model mod'
            elif "D1" in key:
                color = TUBAF.color(ps)['r']
                label = 'model $D_1$'
            elif "A" in key:
                color = TUBAF.color(ps)['o']
                label = 'model $A$'
            elif "HS" in key:
                color = TUBAF.color(ps)['b']
                label = 'HoSi\u2082'
            
            ax.plot(energy[key] - fit_para[key].popt["dE"], fit[key], 
                 label=label, lw=TUBAF.width(ps), color=color)
    ax.plot(en_exp, exp_norm, label='Experiment', color='black', marker='.')

my_plot(ax1, fit)
    
pl.ylim([-0.05,1.95])
pl.xlim([24310,24549])

# pl.legend(bbox_to_anchor=(1., .82),
           # ncol=1, prop={'size':12}, handlelength=1.5)
pl.legend(ncol=2, loc=4, prop={'size':14}, handlelength=1.5, columnspacing=1.)

pl.xlabel('Energy (eV)', fontsize=16)
pl.ylabel('Intensity (a. u.)', fontsize=16)





# ############## inset preedge
# ax2 = pl.axes([0.115,0.475,0.18,0.3], axisbg='white')

# plot_markers(ax2)

# my_plot(ax2, fit)

# pl.setp(ax2.get_xticklabels(), visible=False)
# pl.setp(ax2.get_yticklabels(), visible=False)

# ax2.xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))
# ax2.yaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))

# ax2.get_xaxis().get_major_formatter().set_useOffset(False)
# ax2.get_yaxis().get_major_formatter().set_useOffset(False)

# # limits2 = [8059.5, 8068.7, -0.01, 0.45]
# limits2 = [8061, 8071.5, -0.01, 0.49]
# pl.xlim([limits2[0],limits2[1]])
# pl.ylim([limits2[2],limits2[3]])

# c = TUBAF.color('TUBAF')['r']
# for axis in ['left', 'bottom', 'right', 'top']:
    # ax2.spines[axis].set_color(c)
    # ax2.spines[axis].set_lw(0.7*TUBAF.width(ps))

# def draw_box(ax, limits):
    # """
    # limits = [x1, y1, x2, y2]
    # """
    # x1, x2, y1, y2 = limits
    # ax.plot([x1,x1], [y1,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    # ax.plot([x2,x2], [y1,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    # ax.plot([x1,x2], [y1,y1], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    # ax.plot([x1,x2], [y2,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')

# draw_box(ax1, limits2)

# ax1.annotate("",
            # xy=(8052, 1.44), xycoords='data',
            # xytext=(limits2[0], limits2[3]), textcoords='data',
            # arrowprops=dict(arrowstyle="->",
                            # connectionstyle="arc3",
                            # facecolor=c,
                            # edgecolor=c,
                            # linewidth=0.7*TUBAF.width(ps)), 
            # )






# ############# inset postedge
# ax3 = pl.axes([0.45,0.475,0.18,0.3])
# my_plot(ax3, fit)

# pl.setp(ax3.get_xticklabels(), visible=False)
# pl.setp(ax3.get_yticklabels(), visible=False)

# ax3.xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))
# ax3.yaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))

# ax3.get_xaxis().get_major_formatter().set_useOffset(False)
# ax3.get_yaxis().get_major_formatter().set_useOffset(False)

# # limits3 = [8079.5, 8090, 0.85, 1.35]
# limits3 = [8082.3, 8092.8, 0.9, 1.4]
# pl.xlim([limits3[0],limits3[1]])
# pl.ylim([limits3[2],limits3[3]])

# c = TUBAF.color('TUBAF')['r']
# for axis in ['left', 'bottom', 'right', 'top']:
    # ax3.spines[axis].set_color(c)
    # ax3.spines[axis].set_lw(0.7*TUBAF.width(ps))

# plot_markers(ax3)

# draw_box(ax1, limits3)

# ax1.annotate("",
            # xy=(8100, 1.44), xycoords='data',
            # xytext=(limits3[1], limits3[3]-0.3), textcoords='data',
            # arrowprops=dict(arrowstyle="->",
                            # connectionstyle="arc3",
                            # facecolor=c,
                            # edgecolor=c,
                            # linewidth=0.7*TUBAF.width(ps)), 
            # )





pl.savefig('xafs-compare-K-HoSi2-' + TUBAF.name(ps) + '.pdf')#, transparent=True)
pl.show()

