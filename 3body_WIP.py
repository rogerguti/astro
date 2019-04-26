# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:34:29 2019

@author: Roger Gutierrez Ramon
"""

# Import general packages and modules
import numpy as np
import math as m
import random as rn

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch, Circle
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import proj3d

from scipy.integrate import solve_ivp, odeint

from astropy import constants as const
import spiceypy as spice


# Import personal packages and modules
from astro_equations.solar_system_variables import *
from astro_equations.support_functions_general import set_axes_equal
from astro_equations.support_functions_astro import kepler_third_law_omega, distance, f_potential, three_body, define_initial_conditions


# Define initial conditions for the ode
Y0 = define_initial_conditions()

# Integration
orbit = odeint(three_body, Y0, time)

# Unpacking of variables from ode
sun_x, sun_y, sun_z, asteroid_x, asteroid_y, asteroid_z, sc_x, sc_y, sc_z, sun_vx, sun_vy, sun_vz, asteroid_vx, asteroid_vy, asteroid_vz, sc_vx, sc_vy, sc_vz = orbit.T


# Plotting results and wrapping up

""" fig = plt.figure(1)
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect('equal')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = earth_radius * np.outer(np.cos(u), np.sin(v))
y = earth_radius * np.outer(np.sin(u), np.sin(v))
z = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))

ax.plot_surface(x, y, z, linewidth=0.0, )

ax.plot(moon_x, moon_y, moon_z, color="r", label="Orbit")
set_axes_equal(ax)
plt.show() """


fig = plt.figure(figsize=(9, 9))
ax = fig.add_subplot(111, projection='3d')
ax.plot(asteroid_x, asteroid_y, asteroid_z)
ax.set_aspect('equal')
plt.title('SpiceyPy Cassini Position Example from Jun 20, 2004 to Dec 1, 2005')
plt.show()


#plt.plot(sun_x, sun_y)
# plt.plot(asteroid_x, asteroid_y)
# plt.plot(sc_x, sc_y)
# plt.show()

# plt.plot(time, moon_vx)
# plt.plot(time, moon_vy)
# plt.plot(time, moon_vz)
# plt.show()



# Unloading the Meta-Kernel
spice.unload("./sun_asteroidK.tm")
