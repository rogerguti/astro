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

    sc_initial_x = asteroid_initial_pos_and_vel[0]-sun_initial_pos_and_vel[0]+866
    sc_initial_y = asteroid_initial_pos_and_vel[1]-sun_initial_pos_and_vel[1]+866
    sc_initial_z = asteroid_initial_pos_and_vel[2]-sun_initial_pos_and_vel[2]+866
    sc_initial_vx = asteroid_initial_pos_and_vel[3]-sun_initial_pos_and_vel[3]+0.15
    sc_initial_vy = asteroid_initial_pos_and_vel[4]-sun_initial_pos_and_vel[4]
    sc_initial_vz = asteroid_initial_pos_and_vel[5]-sun_initial_pos_and_vel[5]

    sc_initial_pos_and_vel = [sc_initial_x, sc_initial_y, sc_initial_z, sc_initial_vx, sc_initial_vy, sc_initial_vz]

    Y0 = [0, 0, 0,          # sun position
          asteroid_initial_pos_and_vel[0]-sun_initial_pos_and_vel[0], asteroid_initial_pos_and_vel[1]-sun_initial_pos_and_vel[1], asteroid_initial_pos_and_vel[2]-sun_initial_pos_and_vel[2],             # moon position
          sc_initial_pos_and_vel[0], sc_initial_pos_and_vel[1], sc_initial_pos_and_vel[2],                   # sc position
          0, 0, 0,       # sun velocity
          asteroid_initial_pos_and_vel[3]-sun_initial_pos_and_vel[3], asteroid_initial_pos_and_vel[4]-sun_initial_pos_and_vel[4], asteroid_initial_pos_and_vel[5]-sun_initial_pos_and_vel[5],          # moon velocity
          sc_initial_pos_and_vel[3], sc_initial_pos_and_vel[4], sc_initial_pos_and_vel[5]]                # sc velocity

    return Y0


def eccentricity(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z):

    pos_vector = np.array([pos_x, pos_y, pos_z])
    vel_vector = np.array([vel_x, vel_y, vel_z])

    mag_pos = np.linalg.norm(pos_vector)
    mag_vel = np.linalg.norm(vel_vector)

    evec_num1 = (mag_vel**2 - (mu_sun / mag_pos)) * pos_vector
    evec_num2 = np.dot(pos_vector,vel_vector) * vel_vector

    evec =  (evec_num1 - evec_num2) / mu_sun

    return evec

def lagrangePoints(mu, OP1, OP2):

    # L1 calculation
    coefL1 = [1, 4*mu - 2, 6*mu**2 - 6*mu +1, 4*mu**3 - 6*mu**2 + 4*mu - 1, mu**4 - 2*mu**3 + 5*mu**2 -4*mu + 2, 2*mu**3 - 3*mu**2 + 3*mu - 1]
    L1roots = np.roots(coefL1) 
    L1 = [0,0]
    for i in L1roots:
        if np.isreal(i):
            L1 = [i.real,0]

    # L2 calculation
    coefL2 = [1, 4*mu - 2, 6*mu**2 - 6*mu +1, 4*mu**3 - 6*mu**2 + 2*mu - 1, mu**4 - 2*mu**3 + mu**2 -4*mu + 2, 3*mu - 3*mu**2 - 1]
    L2roots = np.roots(coefL2) 
    L2 = [0,0]
    for i in L2roots:
        if np.isreal(i):
            L2 = [i.real,0]

    # L3 calculation
    coefL3 = [1, 4*mu - 2, 6*mu**2 - 6*mu +1, 4*mu**3 - 6*mu**2 + 2*mu + 1, mu**4 - 2*mu**3 + mu**2 +4*mu - 2, 3*mu**2 - 3*mu + 1]
    L3roots = np.roots(coefL3) 
    L3 = [0,0]
    for i in L3roots:
        if np.isreal(i):
            L3 = [i.real,0]

    # L4 calculation
    L4 = [1/2 - mu, 1/2*m.sqrt(3)]

    # L5 calculation
    L5 = [1/2 - mu, - 1/2*m.sqrt(3)]

    return L1, L2, L3, L4, L5
  
def jacobiCalc(x, y, OP1, OP2):

    r1 = np.sqrt((x+OP1)**2+y**2)
    r2 = np.sqrt((OP2-x)**2+y**2)

    return x**2 + y**2 + 2*OP2/r1 + 2*OP1/r2

def plot3DPotential(mu, OP1, OP2):
    
    LP =  lagrangePoints(mu, OP1, OP2)

    # plot ZVC that run through L1, L2, L3
    delta = 0.025
    x = np.arange(-2.0, 2.0, delta)
    y = np.arange(-2.0, 2.0, delta)
    X, Y = np.meshgrid(x, y)
    Z = -jacobiCalc(X, Y, OP1, OP2)

    for i in range(len(Z)):
        for j in range(len(Z[i])):
            if Z[i][j] < -7:
                Z[i][j] = -7

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax = plt.gca(projection='3d')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    ax.set_aspect('equal')
    #plt.axis('off')
    surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False) # contour for CjL
    #axes.zaxis.set_major_locator(LinearLocator(10))
    #axes.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)

    z_pos_surface = -2.5

    # massive object
    ax.scatter(-OP1, 0, z_pos_surface, s=400, color='r', marker='o', alpha=1)
    #plt.text(-OP1-0.1, -0.15, r'$m_1$')
 
    ax.scatter(OP2, 0, z_pos_surface, s=100, color='g', marker='o', alpha=1)
    #plt.text(OP2-0.1, -0.15, r'$m_2$')
 
    ax.scatter(LP[0][0], LP[0][1], z_pos_surface, s= 50, marker='.', color='b', alpha=0.8)
    ax.scatter(LP[1][0], LP[1][1], z_pos_surface, s= 50, marker='.', color='b', alpha=0.8)
    ax.scatter(LP[2][0], LP[2][1], z_pos_surface, s= 50, marker='.', color='b', alpha=0.8)
    ax.scatter(LP[3][0], LP[3][1], z_pos_surface, s= 50, marker='.', color='b', alpha=0.8)
    ax.scatter(LP[4][0], LP[4][1], z_pos_surface, s= 50, marker='.', color='b', alpha=0.8)
    ax.text(LP[0][0]-0.05, LP[0][1]-0.15, z_pos_surface, r'L$_1$')
    ax.text(LP[1][0]+0.05, LP[1][1]-0.05, z_pos_surface, r'L$_2$')
    ax.text(LP[2][0]-0.2, LP[2][1]-0.05, z_pos_surface, r'L$_3$')
    ax.text(LP[3][0]-0.05, LP[3][1]+0.05, z_pos_surface, r'L$_4$')
    ax.text(LP[4][0]-0.05, LP[4][1]-0.15, z_pos_surface, r'L$_5$')

    plt.title('ZVC for Lagrange points')

    plt.show()

def plotLagrange(mu, OP1, OP2):
    
    LP =  lagrangePoints(mu, OP1, OP2)

    jacobiConst1 = jacobiCalc(LP[0][0], LP[0][1], OP1, OP2)
    jacobiConst2 = jacobiCalc(LP[1][0], LP[1][1], OP1, OP2)
    jacobiConst3 = jacobiCalc(LP[2][0], LP[2][1], OP1, OP2)
    jacobiConst4 = jacobiCalc(LP[3][0], LP[3][1], OP1, OP2)
    jacobiConst5 = jacobiCalc(LP[4][0], LP[4][1], OP1, OP2)

    jacobiConstList = [jacobiConst3, jacobiConst2, jacobiConst1]
    jacobiCilinder = np.sqrt(jacobiConstList)
    colorsString = ['#1E85B8', '#E45B5E', '#E5D7BD']

    # plot ZVC that run through L1, L2, L3
    delta = 0.025
    x = np.arange(-2.0, 2.0, delta)
    y = np.arange(-2.0, 2.0, delta)
    X, Y = np.meshgrid(x, y)
    Z = jacobiCalc(X, Y, OP1, OP2)

    fig, ax = plt.subplots()
    axes = plt.gca()
    axes.set_xlim([-2,2])
    axes.set_ylim([-2,2])
    ax.set_aspect('equal')
    #plt.axis('off')
    plt.contour(X, Y, Z, jacobiConstList, colors=colorsString) # contour for CjL
    jacobiCil1 = Circle((0,0), jacobiCilinder[0], color=colorsString[0], fill=False) 
    jacobiCil2 = Circle((0,0), jacobiCilinder[1], color=colorsString[1], fill=False) 
    jacobiCil3 = Circle((0,0), jacobiCilinder[2], color=colorsString[2], fill=False) 

    ax.add_artist(jacobiCil1)
    ax.add_artist(jacobiCil2)
    ax.add_artist(jacobiCil3)

    # massive object
    CS = plt.scatter(-OP1, 0, marker='o', s=40, color='r', alpha=0.7)
    plt.text(-OP1-0.1, -0.15, r'$m_1$')
 
    CS = plt.scatter(OP2, 0, marker='o', s=20, color='g', alpha=0.7)
    plt.text(OP2-0.1, -0.15, r'$m_2$')
 
    plt.scatter(LP[0][0], LP[0][1], marker='.', color='b', alpha=0.8)
    plt.scatter(LP[1][0], LP[1][1], marker='.', color='b', alpha=0.8)
    plt.scatter(LP[2][0], LP[2][1], marker='.', color='b', alpha=0.8)
    plt.scatter(LP[3][0], LP[3][1], marker='.', color='b', alpha=0.8)
    plt.scatter(LP[4][0], LP[4][1], marker='.', color='b', alpha=0.8)
    CS = plt.text(LP[0][0]-0.05, LP[0][1]-0.15, r'L$_1$')
    CS = plt.text(LP[1][0]+0.05, LP[1][1]-0.05, r'L$_2$')
    CS = plt.text(LP[2][0]-0.2, LP[2][1]-0.05, r'L$_3$')
    CS = plt.text(LP[3][0]-0.05, LP[3][1]+0.05, r'L$_4$')
    CS = plt.text(LP[4][0]-0.05, LP[4][1]-0.15, r'L$_5$')

    plt.title('ZVC for Lagrange points')

    plt.show()

    # plt.savefig('ZVCLagrange.png')  