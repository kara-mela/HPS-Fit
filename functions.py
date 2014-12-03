import os
import evaluationtools as et
from scipy import interpolate
from kara_tools import io
from itertools import product
import numpy as np
import collections

Simulation = collections.namedtuple("Sim", ["E", "Dafs", "Abs"])

def get_dE(keys):
    dE_xafs = et.loaddat('fit-para_enec.dat', todict=True, comment='#')
    dE = {}
    for idx, key in product(dE_xafs, keys):
        if idx.split('_')[0] in key and idx.split('_')[-1] in key:
            dE[key] = dE_xafs[idx][0]
    return dE

def get_exp(R, DIR = os.curdir, norm=None, crop=None):
    if crop==None:
        crop = slice(None,None)
    flist = os.listdir(DIR)
    # fname = filter(lambda s: s.startswith("dafs_hps_%s"%R) and s.endswith("_corr.dat"), flist)
    fname = filter(lambda s: s.startswith("dafs_hps_%s"%R) and s.endswith("_enec.dat"), flist)
    assert len(fname) == 1, "More than 1 file found."
    fname = fname[0]
    data = et.loaddat(fname, todict=True, comment='')
    if norm=="mean":
        data["bragg"] /= data["bragg"][crop].mean()
    elif norm=="max":
        data["bragg"] /= data["bragg"][crop].max()
    return interpolate.interp1d(data["Energy"], data["bragg"], kind="linear")
    
def idx_cut_energy(edge, cut, energy, shift, dE=0):
    """
    determing closest index to cut energy, concerning shift and dE
    """    
    # return list(pl.round_(energy)).index(pl.round_(edge + cut - shift + dE))
    return list(np.round_(energy)).index(np.round_(edge + cut + dE))

def set_k(key_list):
    k = {}
    j = 0
    for key in key_list:
        k[key] = j
        j += 1.
    return k

def patching(fit_F, fit_G, idx_F, idx_G, e_F, e_G):
    """
    ratio for FDM AND Green sim -> adapting to exp
    """
    ratio = fit_F[idx_F] / fit_G[idx_G]
    
    fit_G = np.hstack((fit_F[0:idx_F], fit_G[idx_G:-1]*ratio))
    e_G = np.hstack((e_F[0:idx_F], e_G[idx_G:-1]))
    return fit_G, e_G
    

def get_sim(R, miller, edge):
    """
        loading data from files
        first rough norming to mean()
    """
    Models   = {'A-L23-Green-conv_out_conv.txt'  : 'A_Green', 
                'A-L23-new-all-conv_out_conv.txt': 'A_FDM', 
                'D1-L23-Green-conv_out_conv.txt' : 'D1_Green', 
                'HoSi2-Green-conv_out_conv.txt'  : 'HS_Green', 
                'HoSi2-conv_out_conv.txt'        : 'HS_FDM', 
                'MN-v11_conv.txt'                : 'D1_FDM', 
                'MN-v16_conv.txt'                : 'D1_FDM', 
                # 'modulated-L23-conv_out_conv.txt': 'mod_Green'
                'mod-L3-conv_out.txt': 'mod_Green'}
    result = {}
    for simfile in Models:
        key = "_".join([Models[simfile], R])
        
        #useabs = R=="sat"
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
    return result
