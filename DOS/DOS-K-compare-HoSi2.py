"""
plotting
substance: Ho2PdSi3, Ho2Si4 = HoSi2
edge: L3 (8071eV, X-Ray Data Booklet)
DOS, cols: 4, 12, 24, 40 (s, p, d, f)
"""

from scipy.optimize import curve_fit, leastsq
from scipy.interpolate import interp1d
import pylab as pl
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FixedLocator, MultipleLocator
from kara_tools import TUBAF
from kara_tools import functions as f
import os
# import kara_tools import functions as kf

from matplotlib import rc
rc('font', **{'size':14})

ps = 'TUBA'

edge = 24370 - 4.5
plot_shift = 16 - 4.5# shift of features depending on chosen exp

s, p, d, f, energy = f.get_data(edge)

# Plot fit results
# creating figures, 3 boxes, one for each state, same x-axes
fig = pl.figure()
my_ratios = [1, 1, 1]
gs = gridspec.GridSpec(3,1, height_ratios=my_ratios)
ax = {}
ax["0"] = pl.subplot(gs[0])
ax["1"] = pl.subplot(gs[1], sharex=ax["0"])
ax["2"] = pl.subplot(gs[2], sharex=ax["0"])


# oscillation labels
def plot_markers(ax):
    # feature markers
    my_labels =             [ '$B$', '$C_1$', '$C_2$', '$C_3$', '$C_4$', '$C_5$']
    # my_energies = pl.array( [24.346,  24.360,  24.384,  24.423,  24.480, 24.511])
    my_energies = pl.array( [24.347,  24.363,  24.384,  24.426,  24.483, 24.514])
    my_energies *= 1000
    my_energies += plot_shift
    for i in range(3):
        for line in range(len(my_labels)):  
            if i == 0:
                if line == 0:
                    ax[str(i)].text(my_energies[line]-10, 3.5, my_labels[line], fontsize=16)
                else:
                    ax[str(i)].text(my_energies[line]+.2, 3.5, my_labels[line], fontsize=16)
            
            ax[str(i)].plot([my_energies[line],my_energies[line]], [-1, 20], color='gray', lw=TUBAF.width(ps), linestyle='--')

plot_markers(ax)


for key in s.keys():
    if 'HoSi2' in key:
        mycolor = TUBAF.color(ps)['b']
    elif 'A' in key:
        mycolor = TUBAF.color(ps)['o']
    elif 'D1' in key:
        mycolor = TUBAF.color(ps)['r']
    else:
        mycolor = TUBAF.color(ps)['g']
    
    print key
    # Legende basteln!!!
    if 'HoSi2' in key:
        ax["0"].plot(energy[key], s[key], label='HoSi\u2082', color=mycolor, lw=TUBAF.width(ps))
    elif 'A' in key and '1' in key:
        ax["0"].plot(energy[key], s[key], label='model $A$', color=mycolor, lw=TUBAF.width(ps))
    elif 'D1' in key and '1' in key:
        ax["0"].plot(energy[key], s[key], label='model $D_1$', color=mycolor, lw=TUBAF.width(ps))
    elif 'mod' in key:
    # if 'mod' in key:
        ax["0"].plot(energy[key], s[key], label='model mod', color=mycolor, lw=TUBAF.width(ps))
    else:
        ax["0"].plot(energy[key], s[key], color=mycolor, lw=TUBAF.width(ps))
    
    ax["1"].plot(energy[key], p[key], color=mycolor, lw=TUBAF.width(ps))
    ax["2"].plot(energy[key], d[key], color=mycolor, lw=TUBAF.width(ps))

# scaling
ax["0"].set_ylim(-0.5,4.3)
# ax["0"].yaxis.set_major_locator(MultipleLocator(2))
ax["1"].set_ylim(-0.5,4.3)
# ax["1"].yaxis.set_major_locator(MultipleLocator(2))
ax["2"].set_ylim(-0.5,4.3)
pl.xlim([24310,24549])
ax["0"].xaxis.set_major_locator(FixedLocator((24350, 24400, 24450, 24500)))

bla = ax["0"].get_xticklabels()
for i in bla:
    i.set_fontsize(8)

# ticks and labels
for i in range(3):
    if i != 2:
        ax[str(i)].tick_params(labelbottom='off')

# distance of subplots
pl.subplots_adjust(hspace=0.15)

# ax labels
ax["0"].legend(bbox_to_anchor=(0., 1.1, 1., .065), loc=3,
       ncol=4, mode="expand", borderaxespad=0., prop={'size':14})
       
pl.xlabel('Energy (eV)', fontsize=16)
fig.text(0.04, 0.5, 'LDOS (a. u.)', ha='center', va='center', rotation='vertical', fontsize=16)
ax["0"].set_ylabel('$s$-state', fontsize=14)
ax["1"].set_ylabel('$p$-state', fontsize=14)
ax["2"].set_ylabel('$d$-state', fontsize=14)

# curve labels
"""
the environments of Fig. 3 are NOT Pd-environments!!!
"""
# ax["2"].annotate('(d)', xy=(8057,3.7), xytext=(8043,8.275),arrowprops=dict(arrowstyle="->"), size=12, backgroundcolor='white')
# ax["2"].annotate('(c)', xy=(8056,2.7), xytext=(8043,6.475),arrowprops=dict(arrowstyle="->"), size=12, backgroundcolor='white')
# ax["3"].annotate('(a)', xy=(8055,-.3), xytext=(8043,4.675),arrowprops=dict(arrowstyle="->"), size=12, backgroundcolor='white')
# ax["3"].annotate('(a)', xy=(8056,0.7), xytext=(8043,6.475),arrowprops=dict(arrowstyle="->"), size=12, backgroundcolor='white')
# ax["3"].annotate('(b)', xy=(8057,1.7), xytext=(8043,8.275),arrowprops=dict(arrowstyle="->"), size=12, backgroundcolor='white')

pl.savefig('DOS-K-compare-Hosi2-' + TUBAF.name(ps) + '.pdf', transparent=True)

pl.show()




