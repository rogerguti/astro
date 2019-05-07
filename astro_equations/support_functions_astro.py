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

    # Sun position
    dY[0] = Y[9]
    dY[1] = Y[10]
    dY[2] = Y[11]

    # Asteroid position
    dY[3] = Y[12]
    dY[4] = Y[13]
    dY[5] = Y[14]

    # Spacecraft position
    dY[6] = Y[15]
    dY[7] = Y[16]
    dY[8] = Y[17]

    # distances between objects
    r_sun_asteroid = distance(Y[0], Y[1], Y[2], Y[3], Y[4], Y[5])
    r_sun_sc = distance(Y[0], Y[1], Y[2], Y[6], Y[7], Y[8])
    r_asteroid_sc = distance(Y[3], Y[4], Y[5], Y[6], Y[7], Y[8])

    # Sun velocity
    dY[9] = 0
    dY[10] = 0
    dY[11] = 0
    #dY[9] = f_potential(Y[3], Y[0], r_sun_asteroid, mu_asteroid)
    #dY[10] = f_potential(Y[4], Y[1], r_sun_asteroid, mu_asteroid)
    #dY[11] = f_potential(Y[5], Y[2], r_sun_asteroid, mu_asteroid)

    # Asteroid velocity
    dY[12] = f_potential(Y[0], Y[3], r_sun_asteroid, mu_sun)
    dY[13] = f_potential(Y[1], Y[4], r_sun_asteroid, mu_sun)
    dY[14] = f_potential(Y[2], Y[5], r_sun_asteroid, mu_sun)

    # Spacecraft velocity
    #dY[15] = 0
    #dY[16] = 0
    #dY[17] = 0
    dY[15] = f_potential(Y[0], Y[6], r_sun_sc, mu_sun) + f_potential(Y[3], Y[6], r_asteroid_sc, mu_asteroid)
    dY[16] = f_potential(Y[1], Y[7], r_sun_sc, mu_sun) + f_potential(Y[4], Y[7], r_asteroid_sc, mu_asteroid)
    dY[17] = f_potential(Y[2], Y[8], r_sun_sc, mu_sun) + f_potential(Y[5], Y[8], r_asteroid_sc, mu_asteroid)

    return dY


def define_initial_conditions():
    '''Defines the initial conditions vector

    Input
      Y: state vector for the objects
      t: time vector with the time steps where the ode needs to integrate

    Output
      dY: time derivative of the state vector that needs to be integrated
    '''

    sc_initial_x = asteroid_initial_pos_and_vel[0]+866
    sc_initial_y = asteroid_initial_pos_and_vel[1]+866
    sc_initial_z = asteroid_initial_pos_and_vel[2]+866
    sc_initial_vx = asteroid_initial_pos_and_vel[3]+0.10
    sc_initial_vy = asteroid_initial_pos_and_vel[4]
    sc_initial_vz = asteroid_initial_pos_and_vel[5]

    sc_initial_pos_and_vel = [sc_initial_x, sc_initial_y, sc_initial_z, sc_initial_vx, sc_initial_vy, sc_initial_vz]

    Y0 = [0, 0, 0,          # sun position
          asteroid_initial_pos_and_vel[0], asteroid_initial_pos_and_vel[1], asteroid_initial_pos_and_vel[2],             # moon position
          sc_initial_pos_and_vel[0], sc_initial_pos_and_vel[1], sc_initial_pos_and_vel[2],                   # sc position
          0, 0, 0,       # sun velocity
          asteroid_initial_pos_and_vel[3], asteroid_initial_pos_and_vel[4], asteroid_initial_pos_and_vel[5],          # moon velocity
          sc_initial_pos_and_vel[3], sc_initial_pos_and_vel[4], sc_initial_pos_and_vel[5]]                # sc velocity

    return Y0
