# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 20:58:42 2024

@author: win10
"""

import numpy as np
from sp3file import parse_sp3_file, select_interpolation_data
from sat_pos import sat_pos, clockerror
from ceren_atmos import atmos
from math import sqrt
from numpy.linalg import inv

def calculateReceptionEpoch():
    return (2+2+0+0+6+7+4+0+3+0) * 810

def initializeData():
    approxPosition = np.array([4239146.6414, 2886967.1245, 3778874.4800], dtype=np.float64)
    filePath = r"C:\Users\ceren\Desktop\proje\IGS0OPSFIN_20250910000_01D_15M_ORB.SP3"
    epochs, data = parse_sp3_file(filePath)
    return approxPosition, epochs, data

def getObservations():
    return {
        'G05': 22581252.559,
        'G07': 22830577.076,
        'G08': 24862112.929,
        'G09': 24920406.994,
        'G13': 20928832.115,
        'G14': 20205799.876,
        'G15': 23305396.300,
        'G17': 22705825.247,
        'G19': 24773732.425,
        'G20': 23314693.880,
        'G22': 20423865.035,
        'G30': 21103278.590
    }

def getTGD():
    return {
        'G05': -0.111758708954e-07,
        'G07': -0.111758708954e-07,
        'G08': 0.512227416039e-08,
        'G09': 0.139698386192e-08,
        'G13': -0.111758708954e-07,
        'G14': -0.791624188423e-08,
        'G15': -0.102445483208e-07,
        'G17': -0.111758708954e-07,
        'G19': -0.153668224812e-07,
        'G20': -0.838190317154e-08,
        'G22': -0.111758708954e-07,
        'G30': 0.372529029846e-08
    }

def spp(satelliteObservations, selectedData, approxPosition, epochs, data, receptionEpoch, dayOfYear, alphaParams, betaParams, tgd, includeDelayTGD, c1):
    speedOfLight = 299792458
    receptionEpochWeek = (86400 * 5) + receptionEpoch

    # Calculate satellite positions
    satellitePositions = {}
    for sat, obs in satelliteObservations.items():
        satData = selectedData.get(sat)
        if satData is not None:
            satPos = sat_pos(receptionEpoch, obs, satData, approxPosition)
            satellitePositions[sat] = satPos

    satellitePosMatrix = np.array([satellitePositions[sat] for sat in satellitePositions])

    # Initialize results dictionary
    results = {}
    for i, sat in enumerate(satellitePositions):
        obs = satelliteObservations[sat]
        satData = selectedData.get(sat)
        delay = atmos(dayOfYear, receptionEpoch, receptionEpochWeek, c1, approxPosition, satData, alphaParams, betaParams)
    
        clkData = satData[:, [0, 4]]
        dt = clockerror(receptionEpoch, obs, clkData)
        if includeDelayTGD:
            d = -speedOfLight * dt + delay[4] + delay[5] + delay[3] + (tgd[sat] * speedOfLight)
        else:
            d = -speedOfLight * dt
        results[sat] = {'d': d}

    receiverPosition = approxPosition.copy()
    dValues = np.array([results[sat]['d'] for sat in satellitePositions])
    observationsArray = np.array([satelliteObservations[sat] for sat in satellitePositions])

    # Iterative least squares
    while True:
        A, l = [], []
        for j, sat in enumerate(satellitePositions):
            p0 = np.linalg.norm(satellitePosMatrix[j] - receiverPosition)
            l.append(observationsArray[j] - p0 - dValues[j])
            A.append([(receiverPosition[k] - satellitePosMatrix[j][k]) / p0 for k in range(3)] + [1])

        A = np.array(A)
        l = np.array(l).reshape(-1,1)
        try:
            x = inv(A.T @ A) @ A.T @ l
        except np.linalg.LinAlgError:
            print("Matrix inversion error: singular matrix encountered.")
            break
        dx, dy, dz, dt = x.flatten()
        if all(abs(val) <= 1e-3 for val in (dx, dy, dz)):
            break
        receiverPosition += np.array([dx, dy, dz])

    return receiverPosition

def sppProject():
    receptionEpoch = calculateReceptionEpoch()
    approxPos, epochs, data = initializeData()
    satelliteObservations = getObservations()
    tgd = getTGD()
    c1 = 22719869.219
    dayOfYear = 61
    alphaParams = np.array([0.2794e-07, 0.7451e-08, -0.1192e-06, 0.5960e-07])
    betaParams = np.array([0.1372e+06, -0.3277e+05, -0.6554e+05, -0.5898e+06])
    
    selectedData = select_interpolation_data(dayOfYear, epochs, data)
    
    includeDelayTGD = True  # otomatik evet olarak ayarlandı
    
    updatedPosition = spp(satelliteObservations, selectedData, approxPos, epochs, data, receptionEpoch, dayOfYear, alphaParams, betaParams, tgd, includeDelayTGD, c1)

    print("ESTIMATED COORDINATES")
    print(f"NP: X = {updatedPosition[0]:.3f}, Y = {updatedPosition[1]:.3f}, Z = {updatedPosition[2]:.3f}")

    approxPositionIGS = np.array([4239149.205, 2886968.037, 3778877.204])
    delta = np.abs(updatedPosition - approxPositionIGS)
    msf = np.linalg.norm(updatedPosition - approxPositionIGS)

    print("msf:", msf)
    print("COORDINATES IN IGS RINEX FILE")
    print("*" * 10)
    print(f"X = {approxPositionIGS[0]} m, Y = {approxPositionIGS[1]} m, Z = {approxPositionIGS[2]} m")
    print("*" * 10)
    print("Delta Values:")
    print(f"Delta X = {delta[0]:.3f} m")
    print(f"Delta Y = {delta[1]:.3f} m")
    print(f"Delta Z = {delta[2]:.3f} m")

if __name__ == "__main__":
    sppProject()
