import numpy as np

# speed of light (m/s)

c = 299792458.0  

# earth's angular rotation rate (rad/s)

wE = 7.2921151467e-5  

# approximate position of the receiver in ECEF coordinates (taken from RINEX header)
r_apr = np.array([4239146.6414, 2886967.1245, 3778874.4800])

# Reception epoch: 07:44:00 UTC = 7*3600 + 44*60 = 27840 second
trec = 27840  

# sample pseudorange values (in meters) for two satellites G06 and G11
# these are example values and should be replaced with actual observations if needed

pc_dict = {
    "G06": 22719869.219,  # örnek değer, gerçek değerle değiştir
    "G11": 21509336.237   # örnek değer, gerçek değerle değiştir
}

# precise ephemeris data (SP3 format) for G06, each row represents time(s), X(m), Y(m), Z(m9, clock correction (s)
sp3_G06 = np.array([
    [27000,  6382950.123, 23512300.456, -10876500.789, -0.00025],
    [27100,  6400000.000, 23540000.000, -10850000.000, -0.00023],
    [27200,  6415000.000, 23565000.000, -10830000.000, -0.00021],
    [27300,  6430000.000, 23590000.000, -10810000.000, -0.00020],
    [27400,  6445000.000, 23615000.000, -10790000.000, -0.00019],
    [27500,  6460000.000, 23640000.000, -10770000.000, -0.00018],
    [27600,  6475000.000, 23665000.000, -10750000.000, -0.00017],
    [27700,  6490000.000, 23690000.000, -10730000.000, -0.00016],
    [27800,  6505000.000, 23715000.000, -10710000.000, -0.00015],
    [27900,  6520000.000, 23740000.000, -10690000.000, -0.00014]
])
 # precise ephemeris data (SP3 format) for G11.
sp3_G11 = np.array([
    [27000, 11300000.000, 14380000.000, -17400000.000, -0.00070],
    [27100, 11340000.000, 14395000.000, -17380000.000, -0.00069],
    [27200, 11380000.000, 14410000.000, -17360000.000, -0.00068],
    [27300, 11420000.000, 14425000.000, -17340000.000, -0.00067],
    [27400, 11460000.000, 14440000.000, -17320000.000, -0.00066],
    [27500, 11500000.000, 14455000.000, -17300000.000, -0.00065],
    [27600, 11540000.000, 14470000.000, -17280000.000, -0.00064],
    [27700, 11580000.000, 14485000.000, -17260000.000, -0.00063],
    [27800, 11620000.000, 14500000.000, -17240000.000, -0.00062],
    [27900, 11660000.000, 14515000.000, -17220000.000, -0.00061]
])

# Lagrange interpolation function

def lagrange_interpolation(x, x_vals, y_vals):
    

 
    L = 0.0
    for i in range(len(x_vals)):
        term = y_vals[i]
        for j in range(len(x_vals)):
            if j != i:
                term *= (x - x_vals[j]) / (x_vals[i] - x_vals[j])
        L += term
    return L

# Emission time calculation
"""
    Computes signal emission time using Lagrange interpolation.
    Inputs:
        trec - reception time (seconds)
        pc - pseudorange (meters)
        clk - 10x2 array: [time, clock correction]
    Output:
        tems - emission time (seconds)
 """
def emist(trec, pc, clk):

    t_vals = clk[:, 0]
    clk_vals = clk[:, 1]
    dt_sat = lagrange_interpolation(trec, t_vals, clk_vals)
    tems = trec - (pc / c) - dt_sat
    return tems

# final satellite position calculation

def sat_pos(trec, pc, sp3, r_apr):
    """
   computes the final ECEF coordinates of a satellite at the reception time, 
   correcting for Earth's rotation during signal travel.

   parametres:
   trec: reception time (s)
   pc: pseudorange measurement (m)
   sp3: SP3 satellite data [time, X, Y, Z, clock] (10x5)
   r_apr: approximate receiver position [x, y, z] in meters

   Returns:
   r_sat_final: final satellite coordinates (ECEF, m)
   """
# extarct satellite clock data 

    clk_mat = sp3[:, [0, 4]]
# compute signal emission time

    tems = emist(trec, pc, clk_mat)
# inner function to interpolate satellite coordinates at emission time

    def interpolate_xyz(t_interp, sp3_data):
        pos = []
        for i in range(1, 4):  # X, Y, Z
            pos.append(lagrange_interpolation(t_interp, sp3_data[:, 0], sp3_data[:, i]))
        return np.array(pos)
# approximate satellite coordinates at emission time

    r_sat_apr = interpolate_xyz(tems, sp3)
# compute signal travel time from satellite to receiver

    delta_t = np.linalg.norm(r_sat_apr - r_apr) / c
# compute earths rotation angle durinf that time 

    theta = wE * delta_t
# rotation matrix for correcting earth rotation around z- axis

    R3 = np.array([
        [np.cos(theta), np.sin(theta), 0],
        [-np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])
# apply rotation to get final satellite coordinates at reception time

    r_sat_final = R3 @ r_sat_apr
    return r_sat_final

# compute final satellite coordinates for G06 and G11
final_G06 = sat_pos(trec, pc_dict["G06"], sp3_G06, r_apr)
final_G11 = sat_pos(trec, pc_dict["G11"], sp3_G11, r_apr)

# print the results
print("Final ECEF Coordinates (G06):", final_G06)
print("Final ECEF Coordinates (G11):", final_G11)
