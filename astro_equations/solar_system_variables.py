# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:34:29 2019

@author: Roger Gutierrez Ramon

Groups all variables and constants used for astrodynamics calculations
"""


import numpy as np
import spiceypy as spice


# Load the Meta-Kernel for the planets and asteroid
spice.furnsh("./sun_asteroidK.tm")

# Set the time when the maneouvre will be studied
manoeuvre_time = spice.str2et('2019-10-01T10:00:00')


# Get initial positions and velocities for Sun and asteroid
# All units are meters and seconds
sun_initial_pos_and_vel, sun_ltime = spice.spkezr('SUN', manoeuvre_time, 'ECLIPJ2000',
                                                  'NONE', 'SOLAR_SYSTEM_BARYCENTER')
sun_initial_pos_and_vel *= 1000

asteroid_initial_pos_and_vel, asteroid_ltime = spice.spkezr('2162173', manoeuvre_time, 'ECLIPJ2000',
                                                            'NONE', 'SOLAR_SYSTEM_BARYCENTER')
asteroid_initial_pos_and_vel *= 1000

astronomical_unit = 149597870700
year_in_seconds = 3600 * 24 * 365

integration_time = year_in_seconds / 12
integration_steps = int(integration_time * 2 / 10) # Every 5 seconds

earth_radius = 6371000

m_sun = 1.9885*10**30    # Sun
m_earth = 5.97237*10**24     # Earth
m_moon = 7.342 * 10 ** 22  # Moon
m_sc = 420000  # Spaceship
m_total = m_earth + m_moon

G = 6.672*10**(-11)   # gravitational constant

mu_earth = G * m_earth
mu_moon = G * m_moon
mu_sun = G * m_sun
mu_sc = G * m_sc
mu_asteroid = 30
mu_total = G * m_total

moon_mean_distance = 384399000
moon_mean_tang_velocity = 1018.27
iss_mean_distance = 405000

theta = 0

# Time
time = np.linspace(0, integration_time, integration_steps)
