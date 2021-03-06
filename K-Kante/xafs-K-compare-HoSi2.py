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

ps = 'TUBA'
simplex = False

from matplotlib import rc
rc('font', **{'size':14})

edge = 24365
plot_shift = 16 - 4.5 # shift of features depending on chosen exp

myvars = ["n", "m", "dE"]

data, energy, xafs, fit_para, fit, fitE = {}, {}, {}, {}, {}, {}

# load data, conv 598736
exp_data    = pl.loadtxt('dafs_hps_fluo_pd.dat', skiprows=1)
data['A']   = pl.loadtxt('A-K_conv_out_conv.txt', skiprows=1) 
data['D1']  = pl.loadtxt('D1-K_conv_out_conv.txt', skiprows=1)
data['mod'] = pl.loadtxt('mod-K_out_conv.txt', skiprows=1)

# energy
for key in data:
    energy[key] = data[key][:,0] + edge
    xafs[key] = data[key][:,1]
en_exp = exp_data[:,0]
Exp = interp1d(en_exp, exp_data[:,1], kind='linear') # Energy, flutot

for key in xafs.keys():
    """
    'deleting' sets without intensity
    """
    if xafs[key].max() < 1e-10:
        xafs.pop(key) # No Intensity
    else:
        xafs[key] /= xafs[key].mean()

# fit
print "m, n, dE"
for key in xafs:
    p0 = dict(m=0.0, n=1., c=0., dE=0., Exp=Exp, Isim=xafs[key])
    fit_para[key] = et.fitls(energy[key], pl.zeros(len(energy[key])), 
                         f.Icorr, p0, myvars, fitalg="simplex")
    print fit_para[key].popt["c"]
    
    m, n, dE = fit_para[key].popt["m"], fit_para[key].popt["n"], fit_para[key].popt["dE"]
    print m, n, dE
    
    fit[key] = f.Icorr(energy[key], diff=False, **fit_para[key].popt)
    fitE[key] = energy[key]

R_fact = {}
for key in fit.keys():
    R = key.split('_')[-1]
    # weights = et.gaussian(x=fitE[key], x0=edge, amp=edge, w=1, y0=0.)
    w = 20
    weights = (fitE[key]>(edge-w)) * (fitE[key]<(edge+w))
    R_fact[key] = f.R_factor(exp=Exp(fitE[key]), sim=fit[key], weights=weights)

f.make_fit_dat(fit_para, name='xafs', edge='K', R_fact=R_fact)



# norming
# for key in fit: 
    # fit[key] = (fit[key] - min(fit[key]))/(np.mean(fit[key][:]) - min(fit[key]))
# exp_norm = (Exp(en_exp) - Exp(en_exp[50:]).min())/(np.mean(Exp(en_exp[:450])) - Exp(en_exp[50:]).min())

bla = {'A'  : 800,
       'D1' : 800,
       'mod': -1}
for key in fit: 
    fit[key] = (fit[key] - min(fit[key]))/(np.mean(fit[key][:bla[key]]) - min(fit[key]))
exp_norm = (Exp(en_exp) - Exp(en_exp[50:]).min())/(np.mean(Exp(en_exp[:410])) - Exp(en_exp[50:]).min())
# for key in fit: 
    # fit[key] = (fit[key] - min(fit[key]))/(np.max(fit[key][:700]) - min(fit[key]))
# exp_norm = (Exp(en_exp) - Exp(en_exp[50:]).min())/(np.max(Exp(en_exp[:224])) - Exp(en_exp[50:]).min())




# Plot fit results
# oscillation labels
def plot_markers(ax):
    # feature markers
    my_labels =             ['$B_1$', '$B_2$', '$C_1$', '$C_2$', '$C_3$', '$C_4$']
    # my_energies = pl.array( [24.346,  24.360,  24.384,  24.423,  24.480, 24.511])
    my_energies = pl.array( [ 24.347,  24.363,  24.384,  24.426,  24.483, 24.514])
    my_energies *= 1000
    my_energies += plot_shift
    for i in range(4):
        for line in range(len(my_labels)):  
            if i == 0:
                pl.text(my_energies[line]+.8, 1.85, my_labels[line], fontsize=16)
            
            pl.plot([my_energies[line],my_energies[line]], [-1, 20], 
                     color='gray', lw=0.5*TUBAF.width(ps), linestyle='--')

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


ax1 = pl.axes([.1, .1, .8, .8])

plot_markers(ax1)

my_plot(ax1, fit)
    
# pl.ylim([-0.05,1.15])
pl.ylim([-0.05,1.95])
pl.xlim([24310,24549])

# pl.legend(bbox_to_anchor=(1., .82),
           # ncol=1, prop={'size':12}, handlelength=1.5)
pl.legend(ncol=2, loc=4, prop={'size':14})#, handlelength=1.5, columnspacing=1.)

pl.xlabel('Energy (eV)', fontsize=16)
pl.ylabel('Intensity (a. u.)', fontsize=16)





# ############## inset pd environment
import matplotlib.image as mpimg
ax2 = pl.axes([0.61,0.27,0.27,0.32], axisbg='white')
img = mpimg.imread('Pd-environ-cl.png')
ax2.imshow(img)

pl.setp(ax2.get_xticklabels(), visible=False)
pl.setp(ax2.get_yticklabels(), visible=False)

ax2.xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))
ax2.yaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))

ax2.get_xaxis().get_major_formatter().set_useOffset(False)
ax2.get_yaxis().get_major_formatter().set_useOffset(False)

c = TUBAF.color('TUBAF')['r']
for axis in ['left', 'bottom', 'right', 'top']:
    ax2.spines[axis].set_color(c)
    ax2.spines[axis].set_lw(0.7*TUBAF.width(ps))


    
    
    
    
# ############## inset edge region
ax3 = pl.axes([0.34,0.27,0.25,0.32], axisbg='white')
limits3 = [24370, 24408, 1.26, 1.8]
# pl.ylim([1.18,1.78])
# pl.xlim([24365,24410])
pl.xlim([limits3[0],limits3[1]])
pl.ylim([limits3[2],limits3[3]])

def draw_box(ax, limits, ps):
    """
    limits = [x1, y1, x2, y2]
    """
    if ps == 'TUBAF':
        c = 'black'
    else:
        c = TUBAF.color('TUBAF')['r']
    x1, x2, y1, y2 = limits
    ax.plot([x1,x1], [y1,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    ax.plot([x2,x2], [y1,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    ax.plot([x1,x2], [y1,y1], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    ax.plot([x1,x2], [y2,y2], color=c, lw=0.7*TUBAF.width(ps), linestyle='-')
    
draw_box(ax1, limits3, ps)

my_plot(ax3, fit)

my_energies = pl.array( [ 24.347,  24.363,  24.384,  24.426,  24.483, 24.514])
my_energies *= 1000
my_energies += plot_shift
for i in range(4):
    for line in range(len(my_energies)):  
        pl.plot([my_energies[line],my_energies[line]], [-1, 20], 
                 color='gray', lw=0.5*TUBAF.width(ps), linestyle='--')

ax1.annotate("",
            xy=(24424, 1.17), xycoords='data',
            xytext=(limits3[1], 1.42), textcoords='data',
            arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3",
                            facecolor=c,
                            edgecolor=c,
                            linewidth=0.7*TUBAF.width(ps)), 
            )

pl.setp(ax3.get_xticklabels(), visible=False)
pl.setp(ax3.get_yticklabels(), visible=False)

ax3.xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))
ax3.yaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))

ax3.get_xaxis().get_major_formatter().set_useOffset(False)
ax3.get_yaxis().get_major_formatter().set_useOffset(False)

c = TUBAF.color('TUBAF')['r']
for axis in ['left', 'bottom', 'right', 'top']:
    ax3.spines[axis].set_color(c)
    ax3.spines[axis].set_lw(0.7*TUBAF.width(ps))





pl.savefig('xafs-compare-K-HoSi2-' + TUBAF.name(ps) + '.pdf')#, transparent=True)
pl.show()

