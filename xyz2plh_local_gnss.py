# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 15:02:29 2025

@author: ceren

# 2200674062
# Z.Ceren Vermez
# GNNS Assigment 1

"""
"""
 This code is designed to convert Cartesian coordinates (X, Y, Z) into ellipsoidal coordinates 
 (latitude, longitude, height). The "xyz2plh" function uses an iterative method to calculate latitude.
 The iterations keep going until the difference between two consecutive latitude values is smaller 
 than 10^-12 degrees. This is needed because we don’t know how many iterations it will take to 
 get a precise result. The function updates the latitude and related values in each iteration until 
 the precision is reached.

  Two different methods are used to calculate the zenith angle:
 - The 'atan2' method is more stable for general cases and is preferred for most situations.
 - The 'asin' method gives better results when the zenith angle is close to 90° but can be less stable 
   in some cases. Both methods are used and the results are compared.

 Input:
   x, y, z [Cartesian coordinates]
 Output:
  lat (latitude), lon (longitude), h (height) [Ellipsoidal coordinates]
"""

import math
from math import atan2, degrees, radians, sqrt, sin, cos, asin

# WGS84 Constants
a = 6378137.0          # Semi-major axis [m], WGS84 ellipsoid's semi-major axis
f = 1 / 298.257223563  # Flattening, WGS84 ellipsoid's flattening factor
e2 = 2*f - f**2        # Squared eccentricity, eccentricity of the ellipsoid

def xyz2plh(cart):
    
    """
    Converts Earth-centered (geocentric) Cartesian coordinates (X, Y, Z) to geodetic (ellipsoid) coordinates.
    This conversion uses an iterative approach. The iterations continue until the change in the latitude (𝜑)
    between two successive steps is smaller than 10^-12 degrees.
    
    Input:
        cart: [X, Y, Z] -> Earth-centered Cartesian coordinates (meters)
        
    Output:
        ellp: [lat, lon, h] -> Geodetic coordinates (latitude, longitude, height)
        Latitude (𝜑) and longitude (𝜆) are in degrees, and height (h) is in meters.
    """
    
    x, y, z = cart  # Earth-centered Cartesian coordinates (X, Y, Z)
    lon = atan2(y, x)  # Longitude (𝜆) calculation
    xy = sqrt(x**2 + y**2)  # The distance in the X-Y plane (horizontal distance)
    
    # Initial guess for the latitude (𝜑)
    
    lat = atan2(z, xy * (1 - e2))
    
    # Iterative solution to find the correct latitude
    
    while True:
        N = a / sqrt(1 - e2 * sin(lat)**2)  # Radius of curvature in the prime vertical (N)
        h = xy / cos(lat) - N  # Height (h) calculation
        latitude_new = atan2(z * (N + h), xy * (N * (1 - e2) + h))  # New latitude (𝜑)
        # Stop iterating if the change in latitude is smaller than 10^-12 degrees
        if abs(latitude_new - lat) < 1e-12:
            break
        lat = latitude_new  # Update latitude for the next iteration
    
    # Convert the result from radians to degrees and return
    
    return [degrees(lat), degrees(lon), h]

def local(rec, sat, zenith_method='atan2'):
    
    """
    Computes the azimuth, zenith angles, and slant distance in the local ellipsoidal coordinate system.
    The 'xyz2plh' function is called within the 'local' function to convert the receiver (topocenter) coordinates 
    to geodetic coordinates.
    
    Parameters:
        rec: The receiver (topocenter) coordinates [X, Y, Z] (in meters)
        sat: The satellite (or target) coordinates [X, Y, Z] (in meters)
        zenith_method: 'atan2' (default) or 'asin', method to calculate zenith angle
    
    Output:
        az: Azimuth angle of the target in degrees [0, 360]
        zen: Zenith angle of the target in degrees [-90, 90]
        slantd: Slant range (radial distance) from the topocenter to the target in meters
    """
    
    lat, lon, _ = xyz2plh(rec)  # Convert the receiver's  coordinates to geodetic  coordinates
    lat, lon = radians(lat), radians(lon)  # Convert latitude and longitude from degrees to radians
    
    # Calculate the differences between the satellite and receiver coordinates (dx, dy, dz)
    
    dx, dy, dz = [sat[i] - rec[i] for i in range(3)]
    
    # Transformation matrix components for local coordinate system
    
    slat, clat = sin(lat), cos(lat)  # Calculate sine and cosine of the latitude
    slon, clon = sin(lon), cos(lon)  # Calculate sine and cosine of the longitude
    
    east  = -slat*clon*dx - slon*dy + clat*clon*dz  # Eastward distance component
    north = -slat*slon*dx + clon*dy + clat*slon*dz  # Northward distance component
    up    =  clat*dx + slat*dz  # Upward distance component
    
    # Calculate the azimuth (angle between north and east components)
    
    azimuth = (degrees(atan2(north, east)) + 360) % 360  # Azimuth should be between 0 and 360 degrees
    
    # Calculate the slant distance 
    
    distance = sqrt(dx**2 + dy**2 + dz**2)  # The Euclidean distance between the receiver and the satellite
    
    # Zenith angle calculation 
    
    if zenith_method == 'asin':
        # Method from Leick (2015) GPS Surveying 4th ed.
        zenit = degrees(asin(up / distance))  # Zenith angle in the range of [-90, 90] degrees
    else:
        # Method from Hofmann-Wellenhof (2008) GNSS
        zenit = degrees(atan2(sqrt(east**2 + north**2), up))  # Numerically stable calculation
    
    # Return the computed values: azimuth, zenith angle, and slant range
    
    return azimuth, zenit, distance

# Example Usage: Compute azimuth, zenith angle, and slant distance between Ankara receiver and GPS satellite PRN 5
if __name__ == "__main__":
    # Ankara coordinates (approximate) in ECEF [X, Y, Z] meters
    # Latitude: 39.925533° N, Longitude: 32.866287° E, Height: 850m
    # Converted using https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm
    ankara_receiver = [4207377.0, 2334553.0, 4177103.0]  
    
    # GPS PRN 5 satellite position in ECEF (example from real ephemeris data)
    gps_prn5_satellite = [20200890.5, -13492708.2, 17164787.6]  # [X, Y, Z] meters
    
    # Compute angles and distance
    azimuth, zenith, distance = local(ankara_receiver, gps_prn5_satellite)
    
    # Print the results
    print(f"Ankara (39.9255°N, 32.8663°E) - GPS PRN 5 Satellite:")
    print(f"Azimuth: {azimuth:.4f}°")
    print(f"Zenith Angle: {zenith:.4f}°")
    print(f"Slant Distance: {distance/1000:.3f} km")
