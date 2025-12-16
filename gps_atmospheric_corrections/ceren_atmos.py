from math import radians, degrees
import numpy as np
from sat_pos import emist, sat_pos
from blhandlocal import global2local, xyz2blh
from Ion_Klobuchar import Ion_Klobuchar
from trop_SPP import trop_SPP


trec = 27840  # 7:44
doy = 91 # day of the year
trecw = 86400 *2  + (trec) # GPS week, Sunday+Monday+(epoch)
c1 = 22719869.219 #C1 CODE

rec = np.array([4239146.6414, 2886967.1245, 3778874.4800]) # APPR RECEIVER COORD

sp3 = np.array([
    [6.30,  4238.628861, 23538.892841, 11469.524160, -279.548418],
    [6.45,  3379.513400, 22373.368474, 13864.238163, -279.568227],
    [7.00,  2312.577980, 21044.159569, 16019.736290, -279.588010],
    [7.15,  1033.185865, 19596.634966, 17899.206272, -279.607834],
    [7.30,  -454.541220, 18077.669083, 19470.693115, -279.627624],
    [7.45, -2137.711296, 16533.867238, 20707.621260, -279.647363],
    [8.00, -3994.991451, 15009.822599, 21589.218783, -279.667163],
    [8.15, -5997.263777, 13546.471169, 22100.838747, -279.686954],
    [8.30, -8108.581134, 12179.606466, 22234.174779, -279.706710],
    [8.45, -10287.385102,10938.608750, 21987.369722, -279.726526]
])

alpha = np.array([0.3353e-07, 0.7451e-08, -0.1788e-06, 0.0])

beta = np.array([0.1372e+06, 0.0, -0.3277e+06, 0.2621e+06])


# Function to calculate azimuth, zenith angle, slant distance, ionospheric and tropospheric delays
def atmos(doy, trec, trecw, c1, rec, sp3, alpha, beta):

    # Compute final satellite position
    fpos = sat_pos(trec, c1, sp3, rec)
    
    # Receiver coordinates
    x = rec[0]
    y = rec[1]
    z = rec[2]
    
    # Convert receiver coordinates from XYZ to BLH (latitude, longitude, height)
    blh = xyz2blh(x, y, z)
    
    # Define latitude and longitude in radians
    lon = radians(blh[0])
    lat = radians(blh[1])
    
    # Convert global satellite position to local coordinates
    result = global2local(fpos, rec)
    
    # Azimuth, zenith, and elevation angles in radians
    azm = radians(result[0])
    zen = radians(result[1])
    elv = radians(90) - zen
    
    # Calculate ionospheric delay
    IonD = Ion_Klobuchar(lat, lon, elv, azm, alpha, beta, trecw)
    
    # Use latitude and ellipsoidal height for tropospheric delay calculation
    l1 = blh[1]
    H = blh[2]
    TrD, TrW, ME = trop_SPP(l1, doy, H, elv)
    TrD = ME * TrD
    TrW = ME* TrW
    
    # Convert azimuth and zenith angles back to degrees
    az = result[0]
    zen = result[1]
    slantd = result[2] / 1000 #km
    
    return az, zen, slantd, IonD, TrD, TrW

# Calculate delays
delay = atmos(doy, trec, trecw, c1, rec, sp3, alpha, beta)


print(f"Azimuth  {delay[0]}")
print(f"Zenith {delay[1]}")
print(f"Slant  {delay[2]}")
print(f"Ionospheric delay {delay[3]}")
print(f"Tropospheric dry delay  {delay[4]}")
print(f"Tropospheric wet delay  {delay[5]}")

    
    
    

    
    
    
    
    
    
    


