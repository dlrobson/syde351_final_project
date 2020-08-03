#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.integrate import odeint

# Thermal capacitance of water (J / K * kg)
# https://www.usgs.gov/special-topic/water-science-school/science/heat-capacity-and-water?qt-science_center_objects=0#qt-science_center_objects
# 0.355 kg == 355 mL
c_4 = 4184
m_water = 0.355
C_4 = c_4 * m_water

# Room temperature (K)
T_room = 22 + 273

""" MUG DIMENSIONS (m) """
diameter_inner = 0.12
diameter_outer = 0.15
thickness_stainless_steel = 0.002
length_inner = 0.17
length_outer = 0.20

diameter_air = diameter_outer - 2 * thickness_stainless_steel
thickness_air = diameter_outer - diameter_inner - thickness_stainless_steel
length_air = length_outer - 2 * thickness_stainless_steel

""" INNER STAINLESS STEEL LAYER THERMAL RESISTANCE"""
# Thermal Conductivity of stainless steel type 304 (W / (m * K)), k_stainless_steel
# https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
k_stainless_steel = 14.4
R_5_1 = (
    2
    * thickness_stainless_steel
    / (
        math.pi
        * diameter_inner
        * k_stainless_steel
        * (diameter_inner + 2 * length_inner)
    )
)

""" AIR GAP THERMAL RESISTANCE"""
# Thermal Conductivity of air (W / (m * K)), k_air
# https://www.engineeringtoolbox.com/thermal-conductivity-d_429.html
k_air = 0.0262
R_5_2 = (
    2
    * thickness_air
    / (math.pi * diameter_air * k_air * (diameter_air + 2 * length_air))
)


""" OUTER STAINLESS STEEL LAYER THERMAL RESISTANCE"""
R_5_3 = (
    2
    * thickness_stainless_steel
    / (
        math.pi
        * diameter_outer
        * k_stainless_steel
        * (diameter_outer + 2 * length_outer)
    )
)
# Equivalent resistance
R_5 = R_5_1 + R_5_2 + R_5_3

# Electrical resistance of wire (Ohms), function of temperature
# R_ref - Resistance (Ohms) at 20C
R_ref = 100

# alpha - temp coefficient of resistance for conductor material (Nichrome)
# T_ref - Temperature alpha is specified at (K)
# http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/rstiv.html
alpha = 0.0004
T_ref = 20 + 273

# Voltage of the thermos electrical circuit (V)
V_1 = 9

# function that returns dtemp/dt
def model(temp, t):

    # https://www.cirris.com/learning-center/general-testing/special-topics/177-temperature-coefficient-of-copper
    R_2 = R_ref * (1 + alpha * (temp - T_ref))
    dtemp_dt = -(temp - T_room) / (C_4 * R_5) + (V_1 ** 2) / (C_4 * R_2)
    return dtemp_dt


# IC temp0 in K
temp0 = 100 + 273

# Time points
t = np.linspace(0, 20000)

# Solve ODE
temp = odeint(model, temp0, t)

""" PLOT RESULTS """
# Removes 273 K offset, converting it to C
plt.plot(t, tuple(np.subtract(temp, tuple([273] for _ in range(len(temp))))))
plt.title("Thermos Temperature over time")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (\xb0C)")
plt.show()
