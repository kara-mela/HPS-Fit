import os
import evaluationtools as et
from scipy import interpolate
from kara_tools import io
from itertools import product
import numpy as np
import collections

Simulation = collections.namedtuple("Sim", ["E", "Dafs", "Abs"])

def get_dE(keys):
    dE_xafs = et.loaddat('fit-para-xafs.dat', todict=True, comment='#')
    dE = {}
    for idx, key in product(dE_xafs, keys):
        if idx.split('_')[0] in key and idx.split('_')[-1] in key:
            dE[key] = dE_xafs[idx][1]
    return dE
    
    
def Icorr(E, Isim, Exp, dE=0, m=0, n=1, c=0., Abs=1, diff=True):
    return (m*(E-E[0]) + n) * Isim / Abs**c  - Exp(E - dE) * diff

def get_exp(R, DIR = os.curdir, end="_enec.dat", norm=None, crop=None):
    if crop==None:
        crop = slice(None,None)
    flist = os.listdir(DIR)
    # fname = filter(lambda s: s.startswith("dafs_hps_%s"%R) and s.endswith("_corr.dat"), flist)
    fname = filter(lambda s: s.startswith("dafs_hps_%s"%R) and s.endswith(end), flist)
    assert len(fname) == 1, "More than 1 file found."
    fname = fname[0]
    data = et.loaddat(fname, todict=True, comment='')
    if norm=="mean":
        data["bragg"] /= data["bragg"][crop].mean()
    elif norm=="max":
        data["bragg"] /= data["bragg"][crop].max()
    return interpolate.interp1d(data["Energy"], data["bragg"], kind="linear")

def sort_mod(data):
    """
    for DOS-files
    too many slightly different environments for mod
    -> average
    """
    
    dummy = []
    for key in data.keys():
        if 'mod' not in key:
            continue
        else:
            print key
            dummy.append(data[key])
            data.pop(key, None)
    
    data['mod_Green'] = sum(dummy)/len(dummy)
    
    return data

def sort_D1(data):
    """
    for DOS-files
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
    if edge < 20000: # Ho L edge
        fname  = filter(lambda s: s.startswith("HoSi2-new-out")      and s.endswith("_sd1.txt"), flist) #
        fname += filter(lambda s: s.startswith("HoSi2-Green-out")    and s.endswith("_sd1.txt"), flist) #
        fname += filter(lambda s: s.startswith("A-L23-Green-out_")   and s.endswith("_sd1.txt"), flist) #
        fname += filter(lambda s: s.startswith("A-L23-new-all-out_") and s.endswith("_sd1.txt"), flist) #
        fname += filter(lambda s: s.startswith("D1-L3-sat-out-v2_")  and s.endswith("_sd1.txt"), flist) #
        fname += filter(lambda s: s.startswith("D1-L23-Green-out_")  and s.endswith("_sd1.txt"), flist) #
        fname += filter(lambda s: s.startswith("modulated-L23-out_") and s.endswith("_sd1.txt"), flist) #
    else: # Pd K edge
        fname  = filter(lambda s: s.startswith("A-K-out.txt_")       and s.endswith("_sd0.txt"), flist) # HPS-D-K (extern from Matthias
        fname += filter(lambda s: s.startswith("D1-K-out.txt_")      and s.endswith("_sd0.txt"), flist) # HPS-A-K (extern from Matthias
        fname += filter(lambda s: s.startswith("mod-K_out_")         and s.endswith("_sd0.txt"), flist) # 594415

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
        
        if 'K' in file:
            pass
        elif 'Green' in file:
            key += 'Green_'
        elif 'mod' in file:
            key += 'Green_'
        else:
            key += 'FDM_'
        
        try:
            col = 1 if edge < 20000 else 2
            nr = int(file.split('_')[col])
            nr = str(nr)
        except:
            nr = '1'
        
        key += nr
        
        data[key] = np.loadtxt(file, skiprows=1)
    
    data = sort_mod(data)
    data = sort_D1(data)
    
    for key in data.keys():
        ct = 0.
        
        if 'HoSi2' in key:
            shift = 0.
        elif 'A' in key:
            shift = 1. if edge < 20000 else 0.
            ct = int(key.split('_')[-1]) - 1
        elif 'D1' in key:
            shift = (1. + 2.) if edge < 20000 else 1.
            ct = int(key.split('_')[-1]) - 1
        elif 'mod' in key:
            shift = (1. + 2. + 2.) if edge < 20000 else (1. + 1.)
        
        s[key] = data[key][:,4]  + shift + int(ct)
        p[key] = data[key][:,12] + shift + int(ct)
        d[key] = data[key][:,24] + shift + int(ct)
        if edge < 20000:
            f[key] = data[key][:,40] + shift + int(ct)
        else: f[key] = None
        
        energy[key] = data[key][:,0]
        energy[key] = np.array(energy[key]) + edge
    
    return s, p, d, f, energy
    
def make_fit_dat(fit_para, name='xafs'):
    """
    creating fit-para.dat
    file containing fit parameter err (and dE from xafs)
    """
    header, dE_line, err_line  = [], [], []
    for key in fit_para:
        header.append(key)
        if name == 'xafs':
            dE_line.append(fit_para[key].popt["dE"])
        else:
            dE_line.append(0.0)
        err_line.append(fit_para[key].err) #(residuals(fitted_param)**2).sum()/len(x_m) -> mittlere Fehlerquadrate
    
    content = [err_line, dE_line]
        
    data = dict(zip(header, np.array(content).T))
    et.savedat('fit-para-' + name + '.dat', data, xcol=header[0])
    
def get_sim(R, miller, edge, symbol="L23"):
    """
        loading data from files
        first rough norming to mean()
    """
    if symbol == "L23":
        Models   = {'mod-conv_out.txt'              : 'mod_Green',  # 594937
                    'D1-L23-Green-conv_out_conv.txt': 'D1_Green',   # 580934
                    'D1-L23-conv_conv.txt'          : 'D1_FDM',     # 594935
                    'A-L23-Green-conv_out_conv.txt' : 'A_Green',    # 582128
                    'A-L23-conv_out_conv.txt'       : 'A_FDM',      # 594937
                    'HoSi2-Green-conv_out_conv.txt' : 'HS_Green',   # 580934
                    'HoSi2-conv_out_conv.txt'       : 'HS_FDM'      # 594937
                    }
    elif symbol == "K":
        Models   = {'A-K_conv_out_conv.txt'         : 'A_FDM',      # 595936
                    'D1-K_conv_out_conv.txt'        : 'D1_FDM',     # 594848
                    'mod-K_out_conv.txt'            : 'mod_Green'   # 594848
                    }
        
    result = {}
    for simfile in Models:
        key = "_".join([Models[simfile], R])
        
        # useabs = R=="sat"
        useabs = True
        Ref = miller[R] if not ("HS_" in Models[simfile]) else R
        
        try:
            data = io.FDMNES.loadDAFS(simfile, Ref, absorption=useabs)
        except ValueError: #not found in this model?
            print("Reflection %s not found in file %s"%(R, simfile))
            continue
        
        if useabs:
            Abs = data[2]
        else:
            Abs = np.ones(len(data[0]))
        result[key] = np.vstack((data[0] + edge, 
                                 data[1], 
                                 Abs))
    return result # Energy, Intensity, Abs
