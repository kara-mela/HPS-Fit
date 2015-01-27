"""
data manipulation (adding, norming, ...) and plotting
substance: Ho2PdSi3
edge: L3 (8071eV, X-Ray Data Booklet)
XAFS sim_col: 1
"""

from scipy.optimize import curve_fit, leastsq, fmin
from scipy.interpolate import interp1d
from matplotlib.ticker import FixedLocator, MultipleLocator
import pylab as pl
from itertools import product
import numpy as np
from kara_tools import TUBAF
import evaluationtools as et

ps = 'TUBA'
simplex = False

from matplotlib import rc
rc('font', **{'size':14})

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
plot_shift = 2.7
cut = 42 + plot_shift

data, energy, xafs, fit_para, fit = {}, {}, {}, {}, {}

# load data
data["mod"] = pl.loadtxt("mod-conv-EF125_out.txt", skiprows=1)
data["D1_Green"] = pl.loadtxt("D1-L23-Green-conv_out_conv.txt", skiprows=1)
# data["D1_FDM"] = pl.loadtxt("D1-L23-EF125-conv_conv.txt", skiprows=1)
data["D1_FDM"] = pl.loadtxt("D1-L23-GH40-conv_conv.txt", skiprows=1)
data["A_Green"] = pl.loadtxt("A-L23-Green-conv_out_conv.txt", skiprows=1)
data["A_FDM"] = pl.loadtxt("A-L23-conv-EF125_out_conv.txt", skiprows=1)
data["HS_Green"] = pl.loadtxt("HoSi2-Green-conv_out_conv.txt", skiprows=1)
data["HS_FDM"] = pl.loadtxt("HoSi2-EF125-conv_out_conv.txt", skiprows=1)
exp_data = pl.loadtxt("dafs_hps_average_ho.dat", skiprows=2)

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
    
    dE -= 1.
    fit_para[key][2] -= 1.
    
    fit[key] = Icorr(energy[key], m, n, xafs[key])
    header.append(key)
    if "mod" in key:
        content.append(dE)
    else:
        content.append(dE)
    
data = dict(zip(header, content))
et.savedat('fit-para.dat', data, xcol=header[0])

# norming
# for key in fit: 
    # fit[key] = (fit[key] - min(fit[key]))/(np.median(fit[key]) - min(fit[key]))
# exp_norm = (Exp(en_exp) - Exp(en_exp[23:135]).min())/(Exp(np.median(en_exp[60:135])) - Exp(en_exp[23:135]).min())

for key in fit: 
    fit[key] = (fit[key] - min(fit[key]))/(np.median(fit[key]) - min(fit[key]))
exp_norm = (Exp(en_exp) - Exp(en_exp[140:275]).min())/(Exp(np.median(en_exp[187:275])) - Exp(en_exp[140:275]).min())


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
ax1 = pl.axes([.1, .1, .8, .8])

# oscillation labels
def plot_markers(ax):
    # feature markers
    my_labels =             ['$B_1$', '$B_2$', '$B_3$', '$C_1$', '$C_2$', '$C_3$', '$C_4$']#, '$C_5$', '$C_6$']
    my_energies = pl.array( [ 8044.9, 8061.45,  8065.3,  8074.6,  8084.0,  8103.0,  8138.0])#, 8167, 8192])
    my_energies += plot_shift
    for i in range(4):
        for line in range(len(my_labels)):  
            if i == 0:
                if line == 1 or line == 8 or line == 5:
                    pl.text(my_energies[line]-6, 2.95, my_labels[line], fontsize=16)
                else:
                    pl.text(my_energies[line]+.8, 2.95, my_labels[line], fontsize=16)
            
            pl.plot([my_energies[line],my_energies[line]], [-1, 105], 
                     color='gray', lw=0.5*TUBAF.width(ps), linestyle='--')

plot_markers(ax1)

# border line FDM--Green
# cut_y = 0.02
cut_y = 2.77
pl.plot([edge+cut,edge+cut], [-1, 20], color='.75', lw=TUBAF.width(ps), linestyle='-.')
pl.text(edge+cut+1.5, cut_y+0.03, 'Green', fontsize=14, color='.75')
pl.arrow(edge+cut+2, cut_y, 6, 0., head_width=0.03, head_length=3, fc='.75', ec='.75')
pl.text(edge+cut-9, cut_y+0.03, 'FDM', fontsize=13, color='.75')
pl.arrow(edge+cut-2, cut_y, -4.5, 0., head_width=0.03, head_length=3, fc='.75', ec='.75')

def my_plot(ax, fit):
    for key in fit: 
        if "FDM" not in key:
            print key, ps
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
                label = 'HoSi$_2$'

            ax.plot(energy[key] - fit_para[key][2], fit[key], label=label, lw=TUBAF.width(ps), color=color)
    ax.plot(en_exp, exp_norm, label='Experiment', color='black', marker='.')

my_plot(ax1, fit)
    
pl.ylim([-0.05,3.1])
pl.xlim([8040,8160])
ax1.xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))

# pl.legend(bbox_to_anchor=(1., .82),
           # ncol=1, prop={'size':12}, handlelength=1.5)
pl.legend(ncol=3, loc=4, prop={'size':12}, handlelength=1.5, columnspacing=1.)

pl.xlabel('Energy (eV)', fontsize=16)
pl.ylabel('Intensity (a. u.)', fontsize=16)





############## inset preedge
ax2 = pl.axes([0.115,0.475,0.18,0.3], axisbg='white')

plot_markers(ax2)

my_plot(ax2, fit)

pl.setp(ax2.get_xticklabels(), visible=False)
pl.setp(ax2.get_yticklabels(), visible=False)

ax2.xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))
ax2.yaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))

ax2.get_xaxis().get_major_formatter().set_useOffset(False)
ax2.get_yaxis().get_major_formatter().set_useOffset(False)

# limits2 = [8059.5, 8068.7, -0.01, 0.45]
limits2 = [8061, 8071.5, -0.01, 0.49]
pl.xlim([limits2[0],limits2[1]])
pl.ylim([limits2[2],limits2[3]])

c = TUBAF.color('TUBAF')['r']
for axis in ['left', 'bottom', 'right', 'top']:
    ax2.spines[axis].set_color(c)
    ax2.spines[axis].set_lw(0.7*TUBAF.width(ps))

def draw_box(ax, limits):
    """
    limits = [x1, y1, x2, y2]
    """
    x1, x2, y1, y2 = limits
    ax.plot([x1,x1], [y1,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    ax.plot([x2,x2], [y1,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    ax.plot([x1,x2], [y1,y1], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    ax.plot([x1,x2], [y2,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')

draw_box(ax1, limits2)

ax1.annotate("",
            xy=(8052, 1.44), xycoords='data',
            xytext=(limits2[0], limits2[3]), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3",
                            facecolor=c,
                            edgecolor=c,
                            linewidth=0.7*TUBAF.width(ps)), 
            )






############# inset postedge
ax3 = pl.axes([0.45,0.475,0.18,0.3])
my_plot(ax3, fit)

pl.setp(ax3.get_xticklabels(), visible=False)
pl.setp(ax3.get_yticklabels(), visible=False)

ax3.xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))
ax3.yaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))

ax3.get_xaxis().get_major_formatter().set_useOffset(False)
ax3.get_yaxis().get_major_formatter().set_useOffset(False)

# limits3 = [8079.5, 8090, 0.85, 1.35]
limits3 = [8082.3, 8092.8, 0.9, 1.4]
pl.xlim([limits3[0],limits3[1]])
pl.ylim([limits3[2],limits3[3]])

c = TUBAF.color('TUBAF')['r']
for axis in ['left', 'bottom', 'right', 'top']:
    ax3.spines[axis].set_color(c)
    ax3.spines[axis].set_lw(0.7*TUBAF.width(ps))

plot_markers(ax3)

draw_box(ax1, limits3)

ax1.annotate("",
            xy=(8100, 1.44), xycoords='data',
            xytext=(limits3[1], limits3[3]-0.3), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3",
                            facecolor=c,
                            edgecolor=c,
                            linewidth=0.7*TUBAF.width(ps)), 
            )





pl.savefig('xafs-compare-HoSi2-' + TUBAF.name(ps) + '.pdf')#, transparent=True)
pl.show()

