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
from scipy import integrate as si

ps = 'TUBA' # plotstyle
pl.matplotlib.rc('font', **{'size':14})

density = 7.6062541648443736
edge = 8071 + 6.3
cut = 42
E_lim = slice(0, 350) # only L3 edge
fact=2 # shift of graphs on y axis

# myvars = ["n", "m"] # fit parameters

Reflections = {"sat" : "-215", 
               "110" : "220", 
               "001" : "008",
               "301" : "608"}

ExpFunc = {} # Experimental Data as Functions
Sim = {} # Simulated Data

# k = dict(zip(Reflections.keys(), fact*pl.arange(len(Reflections))))
k = {"sat" : 0., 
     "110" : 1.5, 
     "001" : 3.,
     "301" : 4.5}

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

def Abs_fit(E, Exp, mu, phi, theta, m=0., n=1., d=pl.inf, dE=0., diff=True):
    """
    following Booth2005, equation (1)
    eps_a = 1
    mu_a = 1
    mu_t = mu_f, due to elastical scattering
    """
    g = pl.sin(phi) / pl.sin(theta)
    mf = (m*(E-E[0]) + n) * Isim # linear machine function
    I = mf / (mu + g*mu) * (1 - pl.exp(-abs(d*mu*(1./pl.sin(phi) + 1./pl.sin(theta))))) - Exp(E - dE) * diff
    return I

theta = { # at E = 8075eV
         "sat" : 12.9625957428644,
         "110" : 22.2779615832629, 
         "001" : 11.0648255327577,
         "301" : 43.1643687250282}


# def func(r, mu, r_exp=0):
    # return r**r_exp * pl.exp(-mu*r)
    
# def Ext_func(lambd, mu, r, c_exp=2, br_exp=2, r_exp=0, sin_exp=1):
def Ext_func(lambd, mu, c_exp=2, br_exp=2, r_exp=0, sin_exp=1):
    N = 1.818238e-24
    e = -1.602176565e-19 # C
    m = 9.10938291e-31 # kg
    c = 299792458 # m/s
    # range = (-pl.infty, pl.infty)
    # integral = si.nquad(func=Ext_func, ranges=range, args=(mu=mu, r_exp=r_exp))
    dummy = (N * e**2 / (m * c**c_exp))**br_exp #* lambd**3 / (pl.sin(2*theta)**sin_exp) #* integral
    return dummy

def Ext_fit(E, Exp, mu, phi, theta, g=100., m=0., n=1., d=pl.inf, dE=0., diff=True):
    """
    including correction for extinction following Chandrasekhar1960 http://www.tandfonline.com/doi/pdf/10.1080/00018736000101219
    http://download.springer.com/static/pdf/22/art%253A10.1007%252FBF01596735.pdf?auth66=1424788392_0d0a37d07b6ae8a6ec2c0da883ef423d&ext=.pdf
    see also http://journals.iucr.org/q/issues/1961/11/00/a03314/a03314.pdf
    http://journals.iucr.org/q/issues/1963/11/00/a04006/a04006.pdf
    
    g=0 -> no secondary extinction
    
    to fit: g, phi, theta, m, n (t_0, gamma_0?)
    Laue case: t = t_0/pl.cos(theta) -> I = P_0*Q_0*p_1*r*pl.exp(-mu*r)
    Bragg case: mu*t_0 >> 1:         -> I = P_0*Q_0*p_1/2*mu
        p_1 = (1 + pl.cos(2*theta)**2)/2.
        P_0 = S*I_0 (S is cross section of incident beam)
        Q_0 = abs(N * e**2 * F / (m * c**2)) * lambd**3 / pl.sin(2*theta)
    
    arbitrary crystal shape:
    rho = I_0*Q_0*p_1*V*A(mu_prim)
    mu_prim = mu + g_2*(p_2/p_1)*Q_0
    g_2 = int(W**2) d(theta-theta_B) # theta vs omega!
    W(delta) = distribution function characterizing the misalignment of the mosaic blocks, delta: angular deviation from mean
    p_2 = (1 + pl.cos(2*theta)**4)/2.
    """
    lambd = 12.398 / (E/1000.)
    I_Booth = Abs_fit(E, Exp, mu, phi, theta, m, n, d, dE, diff)
    
    # t_0 = thickness of perfect plate
    # gamma_0 = direction cosine of incident beam relative to normal fo crystal plate
    # -> Zachariasen1945, p. 169
    # g = secondary extinction coefficient = .5 * eta * pl.sqrt(pl.pi)
    # eta = standard deviation of W
    
    # # # alpha     = Ext_func(lambd, mu, r)
    # # # beta_prim = Ext_func(lambd, mu, r, br_exp=4) * (lambd * t_0 / gamma_0)**2 / 3.
    # # # beta_sec  = Ext_func(lambd, mu, r, br_exp=4, c_exp=3, sin_exp=2, r_exp=1) * g * lambd**3
    # # # I = alpha*I_Booth - (beta_prim + beta_sec)*I_Booth**2
    
    # # # # # # # # y = mu*D
    # # # # # # # # x = (N*lambd*F*l)**2
    # # # # # # # # A = pl.exp(-y) * pl.sinh(y) / y # = 1 if already corrected for absorption
    # # # # # # # # B = 1./y - pl.exp(-y) / pl.sinh(y)
    # # # # # # # # if x > 1:
        # # # # # # # # E_L = pl.exp(-y) * pl.sqrt(2./(pl.pi*x)) * (1 - 1./(8*x) - 3./(128*x**2) - 15./(1024*x**3))
    # # # # # # # # else:
        # # # # # # # # E_L = pl.exp(-y) * (1 - x/2. + x**2/4. - 5*x**3/48. + 7*x**4/192.)
    # # # # # # # # E_B = A / pl.sqrt(1 + B*x)
    # # # # # # # # E = E_L*pl.cos(2*theta)**2 + E_B*pl.sin(2*theta)**2
    
    # # # # # # # # I = E * I_Booth
    
    # to-fit: g
    beta = Ext_func(lambd, mu)
    I = I_Booth * (1 - g * beta * I_Booth)
    
    return I


# Fitten
fit_para, fit, fitE = {}, {}, {}
for key in Sim:
    R = key.split("_")[-1]
    E, Isim, Abs = Sim[key]
    
    # p0 = dict(m=0., n=1., theta=theta[R]/180.*pl.pi, phi=theta[R]/180.*pl.pi, d=pl.inf, 
          # Exp=ExpFunc[R], mu=Abs, dE=dE[key])
    # myvars = ["m", "n", "theta", "phi"]
    # fit_para[key] = et.fitls(E, pl.zeros(len(E)), Abs_fit, p0, 
                              # myvars, fitalg='simplex', maxfun=1e6, maxiter=1e6)
    # fit[key] = Abs_fit(E, diff=False, **fit_para[key].popt)
    # fitE[key] = E
    
    p0 = dict(m=0., n=1., theta=theta[R]/180.*pl.pi, phi=theta[R]/180.*pl.pi, d=pl.inf, 
          Exp=ExpFunc[R], mu=Abs, dE=dE[key], g=300.)
    myvars = ["m", "n", "theta", "phi", "g"]
    fit_para[key] = et.fitls(E, pl.zeros(len(E)), Ext_fit, p0, 
                              myvars, fitalg='simplex', maxfun=1e6, maxiter=1e6)
                              
    fit[key] = Ext_fit(E, diff=False, **fit_para[key].popt)
    fitE[key] = E

kf.make_fit_dat(fit_para, name='dafs', edge='L')
        
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
