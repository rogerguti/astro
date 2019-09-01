# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:34:29 2019

@author: Roger Gutierrez Ramon

Groups all variables and constants used for astrodynamics calculations
"""

import numpy as np
import pykep as pk

# General variables
astronomical_unit = 149597870700  # m
G = 6.672 * 10**(-11)  # gravitational constant m3 kg-1 s-2

# Specific variables for space
earth_radius = 6371000
moon_radius = 1737100
sun_radius = 695700000

m_sun = 1.9885 * 10**30  # Sun kg
m_earth = 5.97237 * 10**24  # Earth kg
m_moon = 7.342 * 10**22  # Moon kg
m_sc = 420000  # Spaceship kg

moon_mean_distance = 384399000
moon_mean_tang_velocity = 1018.27
iss_mean_distance = 405000

# Dimensionless variables for 3BP
#mu = 0.2
#mu = 3.04042*10**(-6) # Sun-(Earth+Moon)
mu = 0.0123  # Earth+Moon
#e = 0.054   # Earth+Moon
e = 0.01674 # Sun-(Earth+Moon)

OP1 = mu
OP2 = 1 - mu

# Time
year_in_seconds = 3600 * 24 * 365
month_in_seconds = 3600 * 24 * 30

integration_time = 3 * month_in_seconds

# Initial positions for Sun and EarthMoon System
rot_PrimaryPos = np.array([-OP1, 0, 0])
rot_SecondaryPos = np.array([OP2, 0, 0])

# Variables to adimensionalize other values
distanceBetweenPrimarySecondary = moon_mean_distance
sum_of_masses = m_earth + m_moon
angVel_system = np.sqrt(
    (G * sum_of_masses) / (distanceBetweenPrimarySecondary**3))

# Adimnesionaliziation variables
adimTime = 1 / angVel_system
adimLength = distanceBetweenPrimarySecondary
adimMass = sum_of_masses
adimVel = adimLength / adimTime

adimAngVel = 1

# Variables for the R3BP
primary_radius = earth_radius/adimLength
secondary_radius = moon_radius/adimLength