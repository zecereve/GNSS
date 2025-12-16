# -*- coding: utf-8 -*-
"""
ceren vermez
2200674062


"""

import numpy as np
from lagrange import lagrange
from cal_sp3 import cal_sp3
from rotation import rotation



c = 299792458

w_E = 7.2921151467e-5



def emist(trec, pc, clk):
    
    
    
    t_sat_ce = lagrange(trec, clk) / 1000000
    
    
    
    delta_t = (pc/c)
    tems = trec - delta_t - t_sat_ce
    
    
    return tems



def sat_pos(trec, pc, sp3, r_apr):
    
    
    clk = sp3[:,[0,4]]
    
    tems = emist(trec, pc, clk)
    
    
    
    spos = cal_sp3(tems, sp3)
    
    r_sat = np.array(spos) * 1000
    delta_t = np.linalg.norm(r_sat - r_apr) / c
    
    fpos = rotation(r_sat, w_E*delta_t, 3)
    
    return fpos








    
    
    
    
    
    
    
    
    
    
    