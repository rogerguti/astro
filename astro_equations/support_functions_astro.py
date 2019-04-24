# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:34:29 2019

@author: Roger Gutierrez Ramon

Groups all astro support functions used for astrodynamics calculations
"""


import math as m
from astro_equations.solar_system_variables import *


def kepler_third_law_omega(mu, R):
    '''Calculates the circular angular velocity of an object given the central mass
    parameter and the radius of the orbit

    Input
      mu: central mass paramtere in SI
      R: cicrcular radius of the orbit in SI

    Output
      omega: circular angular velocity in SI
    '''

    omega = m.sqrt(mu/R**3)

    return omega


def distance(ref_x, ref_y, ref_z, target_x, target_y, target_z):
    '''Calculates the distance between target and reference object given the x,y,z cartesian coordinates

    Input
      ref_x: x coordinate of the reference object
      ref_y: y coordinate of the reference object
      ref_z: z coordinate of the reference object
      target_x: x coordinate of the target object
      target_y: y coordinate of the target object
      target_z: z coordinate of the target object

    Output
      d_total: distance between target and reference object
    '''

    d_total = m.sqrt(m.pow(target_x-ref_x, 2) + m.pow(target_y-ref_y, 2) + m.pow(target_z-ref_z, 2))

    return d_total


def f_potential(ref_coord, target_coord, d, mu):
    '''Calculates the gravity potential energy on an object due to another in the direction referenced by the coordinates

    Input
      ref_coord: coordinate of the reference object in the specific axis
      target_coord: coordinate of the target object in the specific axis
      d: total distance between reference and target object
      mu: mass paramtere of the reference object in SI

    Output
      u: gravity potential energy
    '''

    u = -mu * (target_coord - ref_coord) / (m.pow(d, 3))

    return u


def three_body(Y, t):
    '''Main function for the restricted 3 body problem to numerically integrate

    Input
      Y: state vector for the objects
      t: time vector with the time steps where the ode needs to integrate

    Output
      dY: time derivative of the state vector that needs to be integrated
    '''

    dY = [0 for i in range(18)]

    # Earth position
    dY[0] = Y[9]
    dY[1] = Y[10]
    dY[2] = Y[11]

    # Moon position
    dY[3] = Y[12]
    dY[4] = Y[13]
    dY[5] = Y[14]

    # Spacecraft position
    dY[6] = Y[15]
    dY[7] = Y[16]
    dY[8] = Y[17]

    # distances between objects
    r_earth_moon = distance(Y[0], Y[1], Y[2], Y[3], Y[4], Y[5])
    r_earth_sc = distance(Y[0], Y[1], Y[2], Y[6], Y[7], Y[8])
    r_moon_sc = distance(Y[3], Y[4], Y[5], Y[6], Y[7], Y[8])

    # Earth velocity
    dY[9] = f_potential(Y[3], Y[0], r_earth_moon, mu_moon)
    dY[10] = f_potential(Y[4], Y[1], r_earth_moon, mu_moon)
    dY[11] = f_potential(Y[5], Y[2], r_earth_moon, mu_moon)

    # Moon velocity
    dY[12] = f_potential(Y[0], Y[3], r_earth_moon, mu_earth)
    dY[13] = f_potential(Y[1], Y[4], r_earth_moon, mu_earth)
    dY[14] = f_potential(Y[2], Y[5], r_earth_moon, mu_earth)

    # Spacecraft velocity
    dY[15] = f_potential(Y[0], Y[6], r_earth_sc, mu_earth) + f_potential(Y[3], Y[6], r_moon_sc, mu_moon)
    dY[16] = f_potential(Y[1], Y[7], r_earth_sc, mu_earth) + f_potential(Y[4], Y[7], r_moon_sc, mu_moon)
    dY[17] = f_potential(Y[2], Y[8], r_earth_sc, mu_earth) + f_potential(Y[5], Y[8], r_moon_sc, mu_moon)

    return dY


def define_initial_conditions():
    '''Defines the initial conditions vector

    Input
      Y: state vector for the objects
      t: time vector with the time steps where the ode needs to integrate

    Output
      dY: time derivative of the state vector that needs to be integrated
    '''

    center_of_mass_x = moon_mean_distance / (1 + m_earth / m_moon)
    center_of_mass_y = 0
    center_of_mass_z = 0

    mean_ang_velocity = kepler_third_law_omega(mu_total, moon_mean_distance)

    earth_initial_x = - center_of_mass_x
    earth_initial_y = 0
    earth_initial_z = 0
    moon_initial_x = moon_mean_distance - center_of_mass_x
    moon_initial_y = 0
    moon_initial_z = 0
    sc_initial_x = earth_radius - center_of_mass_x + iss_mean_distance
    sc_initial_y = 0
    sc_initial_z = 0
    earth_initial_vx = 0
    earth_initial_vy = mean_ang_velocity * earth_initial_x
    earth_initial_vz = 0
    moon_initial_vx = 0
    moon_initial_vy = mean_ang_velocity * moon_initial_x
    moon_initial_vz = 0
    sc_initial_vx = 0
    sc_initial_vy = moon_mean_tang_velocity*10.5
    sc_initial_vz = 0

    Y0 = [earth_initial_x, earth_initial_y, earth_initial_z,          # earth position
          moon_initial_x, moon_initial_y, moon_initial_z,             # moon position
          sc_initial_x, sc_initial_y, sc_initial_z,                   # sc position
          earth_initial_vx, earth_initial_vy, earth_initial_vz,       # earth velocity
          moon_initial_vx, moon_initial_vy, moon_initial_vz,          # moon velocity
          sc_initial_vx, sc_initial_vy, sc_initial_vz]                # sc velocity

    return Y0
