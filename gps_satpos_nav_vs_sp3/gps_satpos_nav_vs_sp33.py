{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "c4c311aa-3956-4b13-88e7-84a3271cd531",
   "metadata": {},
   "source": [
    "# GMT312 GLOBAL NAVIGATION SATELLITE SYSTEM\n",
    "## Assigment II"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "772b117a-8c72-4eb2-b1d8-fbae97cfae91",
   "metadata": {},
   "source": [
    "### Calculation of the GPS satellite position using the navigation message and precise ephemeris."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "c962607e-12ef-4cea-8859-fde91984813c",
   "metadata": {},
   "source": [
    "### Z. Ceren Vermez\n",
    "### 2200674062"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "3c39da71-ce80-490e-8d6c-d15c1b692e05",
   "metadata": {},
   "source": [
    "### Required Library"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "06048146-3442-4a25-a1c6-d1556eccdbca",
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np\n",
    "from math import sqrt, sin, cos, atan2"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "2464c32e-dfbe-4883-97ed-1c15859f3c92",
   "metadata": {},
   "source": [
    "### 1. First part (Navigation message)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "8c24a738-1deb-4b87-a223-731f3aa27bad",
   "metadata": {},
   "outputs": [],
   "source": [
    "# calculate satellite position from broadcast ephemeris : navigation messages\n",
    "def cal_brd(eph, brd):\n",
    "    \n",
    "    \"\"\"\n",
    "    calculate satellite position from broadcast ephemeris (navigation message)\n",
    "    \n",
    "    Inputs:\n",
    "        eph - calculation epoch in seconds of day (GPS Time)\n",
    "        brd - 7x4 matrix of broadcast orbit parameters\n",
    "        \n",
    "    Output:\n",
    "        spos - 3x1 vector of satellite position (X, Y, Z) in meters\n",
    "    \"\"\"\n",
    "    \n",
    "    sqrt_A = brd[0, 0]        \n",
    "    delta_n = brd[0, 1]       \n",
    "    M0 = brd[0, 2]           \n",
    "    e = brd[0, 3]            \n",
    "    omega = brd[1, 0]         \n",
    "    Cuc = brd[1, 1]           \n",
    "    Cus = brd[1, 2]           \n",
    "    Crc = brd[1, 3]           \n",
    "    Crs = brd[2, 0]           \n",
    "    i0 = brd[2, 1]            \n",
    "    IDOT = brd[2, 2]          \n",
    "    Omega0 = brd[2, 3]        \n",
    "    Omega_dot = brd[3, 0]     \n",
    "    toe = brd[3, 1]           \n",
    "    Cic = brd[4, 0]           \n",
    "    Cis = brd[4, 1]           \n",
    "\n",
    "    mu = 3.986005e14          \n",
    "    omega_e = 7.2921151467e-5 \n",
    "\n",
    "    t = eph\n",
    "    dt = t - toe\n",
    "\n",
    "    A = sqrt_A ** 2\n",
    "    n0 = sqrt(mu / A**3)\n",
    "    n = n0 + delta_n\n",
    "    M = M0 + n * dt\n",
    "\n",
    "    # Kepler's equation solution (simple fixed iteration)\n",
    "    E = M\n",
    "    for _ in range(5):\n",
    "        E = M + e * sin(E)\n",
    "\n",
    "    # true anomaly\n",
    "    v = atan2(sqrt(1 - e**2) * sin(E), cos(E) - e)\n",
    "\n",
    "    # argument of latitude\n",
    "    phi = v + omega\n",
    "\n",
    "    # corrections\n",
    "    du = Cus * sin(2*phi) + Cuc * cos(2*phi)\n",
    "    dr = Crs * sin(2*phi) + Crc * cos(2*phi)\n",
    "    di = Cis * sin(2*phi) + Cic * cos(2*phi)\n",
    "\n",
    "    u = phi + du\n",
    "    r = A * (1 - e * cos(E)) + dr\n",
    "    i = i0 + IDOT * dt + di\n",
    "\n",
    "    x_prime = r * cos(u)\n",
    "    y_prime = r * sin(u)\n",
    "\n",
    "    Omega = Omega0 + (Omega_dot - omega_e) * dt - omega_e * toe\n",
    "\n",
    "    # satellite coordinates\n",
    "    X = x_prime * cos(Omega) - y_prime * cos(i) * sin(Omega)\n",
    "    Y = x_prime * sin(Omega) + y_prime * cos(i) * cos(Omega)\n",
    "    Z = y_prime * sin(i)\n",
    "\n",
    "    spos = np.array([X, Y, Z]).reshape(3, 1)\n",
    "    \n",
    "    return spos\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "d9b1d8d0-27c1-4eb7-b17c-dc957c2952ac",
   "metadata": {},
   "source": [
    "### 2. Second part (Precise ephemeris):"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "4c4032fa-563b-4565-9cb1-343c6b1bae3f",
   "metadata": {},
   "outputs": [],
   "source": [
    "#  lagrange interpolation fonksiyonu\n",
    "def lagrange(eph, dat):\n",
    "    \n",
    "    \"\"\"\n",
    "    9th order Lagrange interpolation.\n",
    "    \n",
    "    Inputs:\n",
    "        eph - interpolation epoch\n",
    "        dat - 10x2 array with time labels (1st column) and variable values (2nd column)\n",
    "        \n",
    "    Output:\n",
    "        out - interpolated value\n",
    "    \"\"\"\n",
    "    \n",
    "    t = dat[:, 0]\n",
    "    y = dat[:, 1]\n",
    "    n = len(t)\n",
    "    out = 0.0\n",
    "\n",
    "    for i in range(n):\n",
    "        term = y[i]\n",
    "        for j in range(n):\n",
    "            if j != i:\n",
    "                term *= (eph - t[j]) / (t[i] - t[j])\n",
    "        out += term\n",
    "    \n",
    "    return out"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "900d9d9f-80b9-4064-a330-631151a1ac13",
   "metadata": {},
   "source": [
    "### 3. Calculate satellite position from precise ephemeris using Lagrange interpolation"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "ed67ccd2-acdd-4f85-a197-d42afad717df",
   "metadata": {},
   "outputs": [],
   "source": [
    "# calculate satellite position feom precise ephemeris usıng 9th order Lagrange ınterpolation.\n",
    "def cal_sp3(eph, sp3):\n",
    "    \n",
    "    \"\"\"\n",
    "    aclculate satellite position from precise ephemeris using Lagrange interpolation.\n",
    "    \n",
    "    Inputs:\n",
    "        eph - calculation epoch in seconds of day\n",
    "        sp3 - 10x4 matrix with time tags and X,Y,Z coordinates\n",
    "        \n",
    "    Output:\n",
    "        spos - 3x1 vector with satellite position (X, Y, Z) in meters\n",
    "    \"\"\"\n",
    "    \n",
    "    t = sp3[:, 0]\n",
    "    X = sp3[:, 1]\n",
    "    Y = sp3[:, 2]\n",
    "    Z = sp3[:, 3]\n",
    "\n",
    "    dat_X = np.column_stack((t, X))\n",
    "    dat_Y = np.column_stack((t, Y))\n",
    "    dat_Z = np.column_stack((t, Z))\n",
    "\n",
    "    X_int = lagrange(eph, dat_X)\n",
    "    Y_int = lagrange(eph, dat_Y)\n",
    "    Z_int = lagrange(eph, dat_Z)\n",
    "\n",
    "    spos = np.array([X_int, Y_int, Z_int]).reshape(3, 1)\n",
    "    \n",
    "    return spos\n"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "8de2a71d-6566-45b2-81d4-fcca4c2ee42b",
   "metadata": {},
   "source": [
    "### Required Data"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "a5671b13-0b90-4092-b611-2b8a4961717d",
   "metadata": {},
   "outputs": [],
   "source": [
    "# example usage with actual data:\n",
    "\n",
    "eph = 23520  # for my school ıd epoch"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "id": "7c16109e-d0d6-44f5-8e55-08779abd8a5e",
   "metadata": {},
   "outputs": [],
   "source": [
    "# broadcast ephemeris data (7x4 matrice)\n",
    "\n",
    "brd_data = np.array([\n",
    "    [82.0, 99.15625, 0.500199406742e-8, 0.788317071176],\n",
    "    [0.516511499882e-5, 0.1055515185e-1, 0.200234353542e-5, 51537.742672],\n",
    "    [180000.0, -0.614672899246e-7, 17.0545511110, 0.141561031342e-6],\n",
    "    [9.48095753479, 332.53125, 4.46210413106, -0.876000774663e-8],\n",
    "    [0.199294015681e-9, 10.0, 23600.0, 0.0],\n",
    "    [20.0, 0.0, 0.465661287308e-8, 82.0],\n",
    "    [179809.0, 40.0, 0.0, 0.0]\n",
    "])\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 20,
   "id": "b66d4836-a55c-4498-93db-00b584a9212e",
   "metadata": {},
   "outputs": [],
   "source": [
    "# precise ephemeris data (10x4 matrice):\n",
    "\n",
    "sp3_data = np.array([\n",
    "    [19800, -26280.254159,  -5125.991542,   1285.187937],\n",
    "    [20700, -26233.567713,  -5377.444618,  -1518.234377],\n",
    "    [21600, -25891.083004,  -5587.628814,  -4296.291687],\n",
    "    [22500, -25254.566484,  -5792.813385,  -7002.638897],\n",
    "    [23400, -24333.918824,  -6028.522771,  -9592.185117],\n",
    "    [24300, -23146.856156,  -6328.041425, -12021.772649],\n",
    "    [25200, -21718.322024,  -6721.006179, -14250.827634],\n",
    "    [26100, -20079.650351,  -7232.134918, -16241.976883],\n",
    "    [27000, -18267.507580,  -7880.135977, -17961.625430],\n",
    "    [27900, -16322.649172,  -8676.837213, -19380.489079]\n",
    "])"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "1a952284-aa5a-4965-abd3-66e7b1a52cfc",
   "metadata": {},
   "source": [
    "### Example Usage:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 22,
   "id": "ec2d418e-1015-40cb-ab75-9a5b7dbf5e5c",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Satellite Position from Broadcast Ephemeris [X, Y, Z] (meters):\n",
      "[[-4080.12372595]\n",
      " [37103.75224738]\n",
      " [-9612.44938835]]\n",
      "\n",
      "Satellite Position from Precise Ephemeris [X, Y, Z] (meters):\n",
      "[[-24190.56269688]\n",
      " [ -6064.08193952]\n",
      " [ -9926.28686594]]\n"
     ]
    }
   ],
   "source": [
    "# calculate position from broadcast ephemeris\n",
    "\n",
    "spos_brd = cal_brd(eph, brd_data)\n",
    "\n",
    "print(\"Satellite Position from Broadcast Ephemeris [X, Y, Z] (meters):\")\n",
    "print(spos_brd)\n",
    "\n",
    "# Calculate position from precise ephemeris\n",
    "\n",
    "spos_sp3 = cal_sp3(eph, sp3_data)\n",
    "\n",
    "print(\"\\nSatellite Position from Precise Ephemeris [X, Y, Z] (meters):\")\n",
    "print(spos_sp3)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "ddfeeb2c-9eaf-46a2-a783-48c8131b7034",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.12.7"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
