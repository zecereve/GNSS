# -*- coding: utf-8 -*-
"""
Created on Sat Apr 26 23:36:51 2025

@author: ceren
"""

# GMT312 GLOBAL NAVIGATION SATELLITE SYSTEM
# Assigment II
# Calculation of the GPS satellite position using the navigation message and precise ephemeris
# Z. Ceren Vermez- 2200674062

#required library

import numpy as np
from math import sqrt, sin, cos, atan2

#1. First part (Navigation message)

# calculate satellite position from broadcast ephemeris : navigation messages
def cal_brd(eph, brd):
    
    """
    calculate satellite position from broadcast ephemeris (navigation message)
    
    Inputs:
        eph - calculation epoch in seconds of day (GPS Time)
        brd - 7x4 matrix of broadcast orbit parameters
        
    Output:
        spos - 3x1 vector of satellite position (X, Y, Z) in meters
    """
    
    sqrt_A = brd[0, 0]        
    delta_n = brd[0, 1]       
    M0 = brd[0, 2]           
    e = brd[0, 3]            
    omega = brd[1, 0]         
    Cuc = brd[1, 1]           
    Cus = brd[1, 2]           
    Crc = brd[1, 3]           
    Crs = brd[2, 0]           
    i0 = brd[2, 1]            
    IDOT = brd[2, 2]          
    Omega0 = brd[2, 3]        
    Omega_dot = brd[3, 0]     
    toe = brd[3, 1]           
    Cic = brd[4, 0]           
    Cis = brd[4, 1]           

    mu = 3.986005e14          
    omega_e = 7.2921151467e-5 

    t = eph
    dt = t - toe

    A = sqrt_A ** 2
    n0 = sqrt(mu / A**3)
    n = n0 + delta_n
    M = M0 + n * dt

    # Kepler's equation solution (simple fixed iteration)
    E = M
    for _ in range(5):
        E = M + e * sin(E)

    # true anomaly
    v = atan2(sqrt(1 - e**2) * sin(E), cos(E) - e)

    # argument of latitude
    phi = v + omega

    # corrections
    du = Cus * sin(2*phi) + Cuc * cos(2*phi)
    dr = Crs * sin(2*phi) + Crc * cos(2*phi)
    di = Cis * sin(2*phi) + Cic * cos(2*phi)

    u = phi + du
    r = A * (1 - e * cos(E)) + dr
    i = i0 + IDOT * dt + di

    x_prime = r * cos(u)
    y_prime = r * sin(u)

    Omega = Omega0 + (Omega_dot - omega_e) * dt - omega_e * toe

    # satellite coordinates
    X = x_prime * cos(Omega) - y_prime * cos(i) * sin(Omega)
    Y = x_prime * sin(Omega) + y_prime * cos(i) * cos(Omega)
    Z = y_prime * sin(i)

    spos = np.array([X, Y, Z]).reshape(3, 1)
    
    return spos

# 2.Second par (precise ephemeris):
    
#  lagrange interpolation fonksiyonu
def lagrange(eph, dat):
    
    """
    9th order Lagrange interpolation.
    
    Inputs:
        eph - interpolation epoch
        dat - 10x2 array with time labels (1st column) and variable values (2nd column)
        
    Output:
        out - interpolated value
    """
    
    t = dat[:, 0]
    y = dat[:, 1]
    n = len(t)
    out = 0.0

    for i in range(n):
        term = y[i]
        for j in range(n):
            if j != i:
                term *= (eph - t[j]) / (t[i] - t[j])
        out += term
    
    return out

# 3. Calculate satellite position from precise ephemeris using Lagrange interpolation

# calculate satellite position feom precise ephemeris usıng 9th order Lagrange ınterpolation.
def cal_sp3(eph, sp3):
    
    """
    aclculate satellite position from precise ephemeris using Lagrange interpolation.
    
    Inputs:
        eph - calculation epoch in seconds of day
        sp3 - 10x4 matrix with time tags and X,Y,Z coordinates
        
    Output:
        spos - 3x1 vector with satellite position (X, Y, Z) in meters
    """
    
    t = sp3[:, 0]
    X = sp3[:, 1]
    Y = sp3[:, 2]
    Z = sp3[:, 3]

    dat_X = np.column_stack((t, X))
    dat_Y = np.column_stack((t, Y))
    dat_Z = np.column_stack((t, Z))

    X_int = lagrange(eph, dat_X)
    Y_int = lagrange(eph, dat_Y)
    Z_int = lagrange(eph, dat_Z)

    spos = np.array([X_int, Y_int, Z_int]).reshape(3, 1)
    
    return spos

# Required datas
# example usage with actual data:

eph = 23520  # for my school ıd epoch

# broadcast ephemeris data (7x4 matrice)

brd_data = np.array([
    [82.0, 99.15625, 0.500199406742e-8, 0.788317071176],
    [0.516511499882e-5, 0.1055515185e-1, 0.200234353542e-5, 51537.742672],
    [180000.0, -0.614672899246e-7, 17.0545511110, 0.141561031342e-6],
    [9.48095753479, 332.53125, 4.46210413106, -0.876000774663e-8],
    [0.199294015681e-9, 10.0, 23600.0, 0.0],
    [20.0, 0.0, 0.465661287308e-8, 82.0],
    [179809.0, 40.0, 0.0, 0.0]
])

# precise ephemeris data (10x4 matrice):

sp3_data = np.array([
    [19800, -26280.254159,  -5125.991542,   1285.187937],
    [20700, -26233.567713,  -5377.444618,  -1518.234377],
    [21600, -25891.083004,  -5587.628814,  -4296.291687],
    [22500, -25254.566484,  -5792.813385,  -7002.638897],
    [23400, -24333.918824,  -6028.522771,  -9592.185117],
    [24300, -23146.856156,  -6328.041425, -12021.772649],
    [25200, -21718.322024,  -6721.006179, -14250.827634],
    [26100, -20079.650351,  -7232.134918, -16241.976883],
    [27000, -18267.507580,  -7880.135977, -17961.625430],
    [27900, -16322.649172,  -8676.837213, -19380.489079]
])

# example usage:
# calculate position from broadcast ephemeris

spos_brd = cal_brd(eph, brd_data)

print("Satellite Position from Broadcast Ephemeris [X, Y, Z] (meters):")
print(spos_brd)

# Calculate position from precise ephemeris

spos_sp3 = cal_sp3(eph, sp3_data)

print("\nSatellite Position from Precise Ephemeris [X, Y, Z] (meters):")
print(spos_sp3)