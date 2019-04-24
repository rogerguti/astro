# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:34:29 2019

@author: Roger Gutierrez Ramon

Groups all variables and constants used for astrodynamics calculations
"""


import numpy as np


astronomical_unit = 149597870700
year_in_seconds = 3600 * 24 * 365

earth_radius = 6371000
m_earth = 5.97237*10**24     # Earth
m_moon = 7.342 * 10 ** 22  # Moon
m_sc = 420000  # Spaceship
m_total = m_earth + m_moon

G = 6.672*10**(-11)   # gravitational constant

mu_earth = G * m_earth
mu_moon = G * m_moon
mu_sc = G * m_sc
mu_total = G * m_total

moon_mean_distance = 384399000
moon_mean_tang_velocity = 1018.27
iss_mean_distance = 405000

theta = 0

# Time
time = np.linspace(0, year_in_seconds, int(year_in_seconds/10))
