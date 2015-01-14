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
import os
# import kara_tools import functions as kf

from matplotlib import rc
rc('font', **{'size':12})

ps = 'TUBAF'

edge = 8071 + 4 
cut  = 45

def sort_mod(data):
    """
    too many slightly different environments for mod
    -> average
    """
    
    dummy = []
    for key in data.keys():
        if 'mod' not in key:
            continue
        else:
            dummy.append(data[key])
            data.pop(key, None)
    
    data['mod_Green'] = sum(dummy)/32.
    
    return data

def sort_D1(data):
    """
    environments forming two groups
    choosing one representative of each
    """
    
    for key in data.keys():
        if 'D1' not in key:
            continue
        else:
            if int(key.split('_')[-1]) > 2:
                data.pop(key, None)
    
    return data

def get_data(edge, DIR = os.curdir):
    """
    search directory DIR for DOS-files
    give keys (Green/FDM, model)
    return DOS of orbitals s, p, d, f
    return energy
    """
    
    flist = os.listdir(DIR)
    fname  = filter(lambda s: s.startswith("HoSi2-new-out")      and s.endswith("_sd1.txt"), flist)
    fname += filter(lambda s: s.startswith("HoSi2-Green-out")    and s.endswith("_sd1.txt"), flist)
    fname += filter(lambda s: s.startswith("A-L23-Green-out_")   and s.endswith("_sd1.txt"), flist)
    fname += filter(lambda s: s.startswith("A-L23-new-all-out_") and s.endswith("_sd1.txt"), flist)
    fname += filter(lambda s: s.startswith("D1-L3-sat-out-v2_")  and s.endswith("_sd1.txt"), flist)
    fname += filter(lambda s: s.startswith("D1-L23-Green-out_")  and s.endswith("_sd1.txt"), flist)
    fname += filter(lambda s: s.startswith("modulated-L23-out_") and s.endswith("_sd1.txt"), flist)

    data = {}
    s, p, d, f = {}, {}, {}, {}
    energy = {}
    
    for file in fname:
        key = ''
        if file.startswith('H'):
            key += 'HoSi2_'
        elif file.startswith('A'):
            key += 'A_'
        elif file.startswith('D'):
            key += 'D1_'
        elif file.startswith('m'):
            key += 'mod_'
        
        if 'Green' in file:
            key += 'Green_'
        elif 'mod' in file:
            key += 'Green_'
        else:
            key += 'FDM_'
        
        key += file.split('_')[1]
        
        data[key] = pl.loadtxt(file, skiprows=1)
    
    data = sort_mod(data)
    data = sort_D1(data)
    
    for key in data.keys():
        ct = 0.
        
        if 'HoSi2' in key:
            shift = 0.
        elif 'A' in key:
            shift = 1.
            ct = int(key.split('_')[-1]) - 1
        elif 'D1' in key:
            shift = 1. + 2.
            ct = int(key.split('_')[-1]) - 1
        elif 'mod' in key:
            shift = 1. + 2. + 2.
        
        s[key] = data[key][:,4]  + shift + int(ct)
        p[key] = data[key][:,12] + shift + int(ct)
        d[key] = data[key][:,24] + shift + int(ct)
        f[key] = data[key][:,40] + shift + int(ct)
        
        energy[key] = data[key][:,0]
        energy[key] = pl.array(energy[key]) + edge
    
    return s, p, d, f, energy

s, p, d, f, energy = get_data(edge)

# Daten beschneiden und zusammenfuehren
for key in s.keys():
    if not "FDM" in key:
        continue
    
    keyG = key.replace("FDM", "Green")
    
    idxF = s[key][0] <= (edge + cut)
    idxG = s[keyG][0] > (edge + cut)
    
    ratio = s[keyG][idxG] / s[key][~idxF]
    s[keyG][idxG] /= ratio
    ratio = p[keyG][idxG] / p[key][~idxF]
    p[keyG][idxG] /= ratio
    ratio = d[keyG][idxG] / d[key][~idxF]
    d[keyG][idxG] /= ratio
    ratio = f[keyG][idxG] / f[key][~idxF]
    f[keyG][idxG] /= ratio
    
    s[key] = pl.hstack((s[key][:~idxF], s[keyG][idxG:]))
    p[key] = pl.hstack((p[key][:~idxF], p[keyG][idxG:]))
    d[key] = pl.hstack((d[key][:~idxF], d[keyG][idxG:]))
    f[key] = pl.hstack((f[key][:~idxF], f[keyG][idxG:]))
    
    s.pop(keyG)
    p.pop(keyG)
    d.pop(keyG)
    f.pop(keyG)

# Plot fit results
# creating figures, 4 boxes, one for each state, same x-axes
fig = pl.figure()
my_ratios = [1, 1, 1, 1]
gs = gridspec.GridSpec(4,1, height_ratios=my_ratios)
ax = {}
ax["0"] = pl.subplot(gs[0])
ax["1"] = pl.subplot(gs[1], sharex=ax["0"])
ax["2"] = pl.subplot(gs[2], sharex=ax["0"])
ax["3"] = pl.subplot(gs[3], sharex=ax["0"])

for key in s.keys():
    if 'HoSi2' in key:
        mycolor = TUBAF.color(ps)['b']
    elif 'A' in key:
        mycolor = TUBAF.color(ps)['o']
    elif 'D1' in key:
        mycolor = TUBAF.color(ps)['r']
    else:
        mycolor = TUBAF.color(ps)['g']
    
    # Legende basteln!!!
    if 'HoSi2' in key:
        ax["0"].plot(energy[key], s[key], label='HoSi$_2$', color=mycolor, lw=TUBAF.width(ps))
    elif 'A' in key and '2' in key:
        ax["0"].plot(energy[key], s[key], label='model $A$', color=mycolor, lw=TUBAF.width(ps))
    elif 'D1' in key and '2' in key:
        ax["0"].plot(energy[key], s[key], label='model $D_1$', color=mycolor, lw=TUBAF.width(ps))
    elif 'mod' in key:
        ax["0"].plot(energy[key], s[key], label='model mod', color=mycolor, lw=TUBAF.width(ps))
    else:
        ax["0"].plot(energy[key], s[key], color=mycolor, lw=TUBAF.width(ps))
    
    ax["0"].plot(energy[key], s[key], color=mycolor, lw=TUBAF.width(ps))
    ax["1"].plot(energy[key], p[key], color=mycolor, lw=TUBAF.width(ps))
    ax["2"].plot(energy[key], d[key], color=mycolor, lw=TUBAF.width(ps))
    ax["3"].plot(energy[key], f[key], color=mycolor, lw=TUBAF.width(ps))

# scaling
ax["0"].set_ylim(-0.5,9.9)
ax["0"].yaxis.set_major_locator(MultipleLocator(2))
ax["1"].set_ylim(-0.5,9.9)
ax["1"].yaxis.set_major_locator(MultipleLocator(2))
ax["2"].set_ylim(-0.5,9.9)
ax["3"].set_ylim(-0.5,9.9)
ax["0"].set_xlim([8040,8160])
ax["0"].xaxis.set_major_locator(FixedLocator((8050, 8075, 8100, 8125, 8150)))

bla = ax["0"].get_xticklabels()
for i in bla:
    i.set_fontsize(8)

# ticks and labels
for i in range(4):
    if i != 3:
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
ax["3"].set_ylabel('$f$-state', fontsize=14)

# curve labels
ax["2"].annotate('(d)', xy=(8057,3.7), xytext=(8043,8.4),arrowprops=dict(arrowstyle="->"), size=12)
ax["2"].annotate('(c)', xy=(8056,2.7), xytext=(8043,6.6),arrowprops=dict(arrowstyle="->"), size=12)
ax["3"].annotate('(a)', xy=(8055,-.3), xytext=(8043,4.8),arrowprops=dict(arrowstyle="->"), size=12, backgroundcolor='white')
ax["3"].annotate('(a)', xy=(8056,0.7), xytext=(8043,6.6),arrowprops=dict(arrowstyle="->"), size=12)
ax["3"].annotate('(b)', xy=(8057,1.7), xytext=(8043,8.4),arrowprops=dict(arrowstyle="->"), size=12)

# base = 3.5
# fact = 0.7
# ax["2"].text(8161, base+3*fact, '(d)', fontsize=11, color=TUBAF.color(ps)['r'])
# ax["2"].text(8161, base+1*fact, '(c)', fontsize=11, color=TUBAF.color(ps)['r'])
# ax["2"].text(8161, base-1*fact, '(b)', fontsize=11, color=TUBAF.color(ps)['o'])
# ax["2"].text(8161, base-3*fact, '(a)', fontsize=11, color=TUBAF.color(ps)['o'])
# ax["2"].text(8161, base-5*fact, '(a)', fontsize=11, color=TUBAF.color(ps)['b'])

# oscillation labels
def plot_markers(ax):
    # feature markers
    my_labels =             ['$B_1$', '$B_2$', '$C_1$', '$C_2$', '$C_3$']#, '$C_4$', '$C_5$']
    my_energies = pl.array( [8061.4, 8065.3, 8074.6, 8103, 8138])#, 8167, 8192])
    for i in range(4):
        for line in range(len(my_labels)):  
            if i == 0:
                if line == 0 or line == 6:
                    ax[str(i)].text(my_energies[line]-6, 7.9, my_labels[line], fontsize=16)
                else:
                    ax[str(i)].text(my_energies[line]+.2, 7.9, my_labels[line], fontsize=16)
            
            ax[str(i)].plot([my_energies[line],my_energies[line]], [-1, 105], color='gray', lw=TUBAF.width(ps), linestyle='--')

plot_markers(ax)


# # # border line FDM--Green
# # pl.plot([edge+cut,edge+cut], [-1, 105], color='gray', lw=TUBAF.width(ps), linestyle='--')
# # pl.text(edge+cut+5.5, 1., 'Green', fontsize=14, color='0.33')
# # pl.arrow(edge+cut+5, .97, 10, 0., head_width=0.02, head_length=5, fc='0.33', ec='0.33')
# # pl.text(edge+cut-13, 1., 'FDM', fontsize=14, color='0.33')
# # pl.arrow(edge+cut-5, .97, -10, 0., head_width=0.02, head_length=5, fc='0.33', ec='0.33')
 

pl.savefig('DOS-compare-Hosi2-' + TUBAF.name(ps) + '.pdf', transparent=True)

pl.show()




