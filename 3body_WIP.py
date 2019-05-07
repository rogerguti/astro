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
from astro_equations.support_functions_astro import kepler_third_law_omega, distance, f_potential, three_body, define_initial_conditions, eccentricity


# Define initial conditions for the ode
Y0 = define_initial_conditions()

print("len time = " , len(time))
print("time = " , time[0], ", " , time[-1])

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


time1 = 0
time2 = int(len(time)/2)
time3 = len(time) - 1

test1 = eccentricity(asteroid_x[time1], asteroid_y[time1], asteroid_z[time1], asteroid_vx[time1], asteroid_vy[time1], asteroid_vz[time1])
test2 = eccentricity(asteroid_x[time2], asteroid_y[time2], asteroid_z[time2], asteroid_vx[time2], asteroid_vy[time2], asteroid_vz[time2])
test3 = eccentricity(asteroid_x[time3], asteroid_y[time3], asteroid_z[time3], asteroid_vx[time3], asteroid_vy[time3], asteroid_vz[time3])

magtest1 = np.linalg.norm(test1)
magtest2 = np.linalg.norm(test2)
magtest3 = np.linalg.norm(test3)

print("magtest1 = " , magtest1)
print("magtest2 = " , magtest2)
print("magtest3 = " , magtest3)
print("len time = " , len(time))
print("time = " , time[0], ", " , time[-1])

'''
r_max = 0
r_min = 100000000000000000000000000000000000000

for i in range(len(asteroid_x)):

    r_new = distance(sun_x[i], sun_y[i], sun_z[i], asteroid_x[i], asteroid_y[i], asteroid_z[i])

    if r_new > r_max:
        r_max = r_new
    if r_new < r_min:
        r_min = r_new


eccentricity = m.sqrt(1-(r_min**2/r_max**2))

print("r_max: ", r_max)
print("r_min: ", r_min)
print("eccentricity: ", eccentricity)

'''

orbit_sc_x = sc_x - asteroid_x
orbit_sc_y = sc_y - asteroid_y
orbit_sc_z = sc_z - asteroid_z

ryugu_circleplot1 = plt.Circle((0, 0), 500, color='r')
ryugu_circleplot2 = plt.Circle((0, 0), 500, color='r')
ryugu_circleplot3 = plt.Circle((0, 0), 500, color='r')

fig1 = plt.figure(1, figsize=(9, 9))
ax = fig1.add_subplot(111, projection='3d')
ax.plot(orbit_sc_x, orbit_sc_y, orbit_sc_z)
# draw sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 500 * np.outer(np.cos(u), np.sin(v))
y = 500 * np.outer(np.sin(u), np.sin(v))
z = 500 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='r')
#ax.plot(sc_x, sc_y, sc_z)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_aspect('equal')
plt.title('Orbit of the marker around Ryugu during ' + str(integration_time/ (3600 * 24)) + " days")

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(131)
ax2.plot(orbit_sc_x, orbit_sc_y)
ax2.add_artist(ryugu_circleplot1)
ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_aspect('equal')
ax2 = fig2.add_subplot(132)
ax2.plot(orbit_sc_x, orbit_sc_z)
ax2.add_artist(ryugu_circleplot2)
ax2.set_xlabel('X')
ax2.set_ylabel('Z')
ax2.set_aspect('equal')
ax2 = fig2.add_subplot(133)
ax2.plot(orbit_sc_y, orbit_sc_z)
ax2.add_artist(ryugu_circleplot3)
ax2.set_xlabel('Y')
ax2.set_ylabel('Z')
ax2.set_aspect('equal')


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
