# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:34:29 2019

@author: Roger Gutierrez Ramon

Groups all astro support functions used for astrodynamics calculations
"""

import numpy as np
import math as m

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
from matplotlib import animation
from matplotlib.animation import FFMpegWriter
from astro_equations.variables_3bp import *


class FasterFFMpegWriter(FFMpegWriter):
    """
      FFMpeg-pipe writer bypassing figure.savefig.
    """
    def __init__(self, **kwargs):
        """Initialize the Writer object and sets the default frame_format."""
        super().__init__(**kwargs)
        self.frame_format = 'argb'

    def grab_frame(self, **savefig_kwargs):
        """Grab the image information from the figure and save as a movie frame.

        Doesn't use savefig to be faster: savefig_kwargs will be ignored.
        """
        try:
            # re-adjust the figure size and dpi in case it has been changed by the
            # user.  We must ensure that every frame is the same size or
            # the movie will not save correctly.
            self.fig.set_size_inches(self._w, self._h)
            self.fig.set_dpi(self.dpi)
            # Draw and save the frame as an argb string to the pipe sink
            self.fig.canvas.draw()
            self._frame_sink().write(self.fig.canvas.tostring_argb())
        except (RuntimeError, IOError) as e:
            out, err = self._proc.communicate()
            raise IOError('Error saving animation to file (cause: {0}) '
                          'Stdout: {1} StdError: {2}. It may help to re-run '
                          'with --verbose-debug.'.format(e, out, err))


def kepler_third_law_ang_vel(mu, r):
    """Calculates the circular angular velocity of the Kepler 2BP.
    
    Calculates the circular angular velocity of an object in the Kepler 2BP 
    given the central mass parameter and the radius of the orbit.

    Args:
      mu: Central mass paramtere in SI.
      r: Cicrcular radius of the orbit in SI.

    Returns:
      omega: Circular angular velocity in SI.
    """

    omega = np.sqrt(mu / r**3)

    return omega


def distance(ref_object, Y, mu):
    """Calculates the distance between target and reference object.

    Args:
      ref_object (str): Which reference object, 'main' or 'secondary'.
      Y: State vector for the spacecraft.
      mu: Mass parameter of the system.
      
    Returns:
      r: Distance between target and reference object.
    """

    if ref_object == 'main':

        r = np.sqrt((mu + Y[0])**2 + Y[1]**2 + Y[2])

    elif ref_object == 'secondary':

        r = np.sqrt((Y[0] - 1 + mu)**2 + Y[1]**2 + Y[2]**2)

    return r


def potential_dU(coord, Y, mu):
    """Calculates the derivative of the scalar potential U in each coordinate.
    
    Calculates the derivative of the scalar potential U that accounts for 
    gravitational and centrifugal forces wrt the variable in 'coord'.

    Args:
      coord: Variable wrt which to derivate the scalar potential U.
      Y: State vector for the spacecraft.
      mu: Mass constant of the system.

    Returns:
      du: Derivative of the potential U.
    """

    if coord == 'x':

        du = Y[0] - (1 - mu) * (Y[0] + mu) / (distance(
            'main', Y, mu)**3) - mu * (Y[0] - 1 + mu) / (distance(
                'secondary', Y, mu)**3)

    elif coord == 'y':

        du = Y[1] - (1 - mu) * Y[1] / (distance(
            'main', Y, mu)**3) - mu * Y[1] / (distance('secondary', Y, mu)**3)

    elif coord == 'z':

        du = Y[2] - (1 - mu) * Y[2] / (distance(
            'main', Y, mu)**3) - mu * Y[2] / (distance('secondary', Y, mu)**3)

    return du


def state_vector_CR3BP(t, Y):
    """ODE for the circular restricted 3 body problem to numerically integrate.

    Args:
      Y: State vector for the objects.
      t: Time vector with the time steps where the ode needs to integrate.

    Returns:
      dY: Time derivative of the state vector that needs to be integrated.
    """

    dY = [0 for i in range(6)]
    
    # Spacecraft position in rotating frame
    dY[0] = Y[3]
    dY[1] = Y[4]
    dY[2] = Y[5]

    # Spacecraft velocity
    dY[3] = 2 * Y[4] + potential_dU('x', Y, mu)
    dY[4] = -2 * Y[3] + potential_dU('y', Y, mu)
    dY[5] = -Y[2] + potential_dU('z', Y, mu)

    return dY

def state_vector_ER3BP(f, Y):
    """ODE for the elliptic restricted 3 body problem to numerically integrate.

    Args:
      Y: State vector for the objects.
      f: True anomaly vector with the true anomaly steps where the ode needs to integrate.

    Returns:
      dY: Derivative respect to the true anomaly of the state vector that needs to be integrated.
    """

    dY = [0 for i in range(6)]

    # Spacecraft position in rotating frame
    dY[0] = Y[3]
    dY[1] = Y[4]
    dY[2] = Y[5]

    # Spacecraft velocity
    dY[3] = 2 * Y[4] + (1 / (e * np.cos(f)) * potential_dU('x', Y, mu)
    dY[4] = -2 * Y[3] + (1 / (e * np.cos(f)) * potential_dU('y', Y, mu)
    dY[5] = -Y[2] + (1 / (e * np.cos(f)) * potential_dU('z', Y, mu)

    return dY


def define_initial_conditions(mu):
    """Defines the initial coordinates, time and other optional graphics bools.

    Args:
      mu: Mass parameter of the CR3BP.

    Returns:
      Y0: Initial conditions vector for the SC, XYZ VX VY VZ.
      timeSpan: Tange of time to integrate.
      timePoints: Times at which to store the computed solution.
      boolZVC (bool): Boolean to plot the ZVC and Lagrange Points.
      boolLagrange (bool): Boolean to plot only the Lagrange Points.
      boolAnimation (list): List with information of the data as animation:

        - (bool): True to create animation, otherwise False.

        - (bool): True to save the animation, otherwise False.

        - Type of animation (str): 'Standard' or 'Comet'-like.

        - Reference System (str): 'Rotating' or 'Inertial'.
    """

    rot_SCPos_initial = np.array([1 - mu, 0.0455,
                                  0.0])  # SC position in rotating frame
    rot_SCVel_initial = np.array([-0.5322, 0.2,
                                  0.0])  # SC velocity in rotating frame
    #rot_SCVel_initial = np.array([-0.4, 0.2, 0.0]) # SC velocity in rotating frame

    Y0 = np.concatenate(
        (rot_SCPos_initial, rot_SCVel_initial),
        axis=None)  # SC position and velocity in rotating frame

    timeSpan = [0, int(integration_time / adimTime)]
    timePoints = np.linspace(timeSpan[0], timeSpan[1], 15000)

    boolLagrange = False
    boolZVC = True
    boolAnimation = [
        False, True, 'Comet', 'Inertial'
    ]  # Animation True/False, Save True/False, Type Standard/Comet, Ref.Frame
    # Rotating/Inertial

    return Y0, timeSpan, timePoints, boolZVC, boolLagrange, boolAnimation


def eccentricity(pos_x, pos_y, pos_z, vel_x, vel_y, vel_z):

    pos_vector = np.array([pos_x, pos_y, pos_z])
    vel_vector = np.array([vel_x, vel_y, vel_z])

    mag_pos = np.linalg.norm(pos_vector)
    mag_vel = np.linalg.norm(vel_vector)

    evec_num1 = (mag_vel**2 - (mu_Sun / mag_pos)) * pos_vector
    evec_num2 = np.dot(pos_vector, vel_vector) * vel_vector

    evec = (evec_num1 - evec_num2) / mu_Sun

    return evec


def calculate_lagrange_points(mu):
    """Calculates the position of the Lagrange Points of the CR3BP in the 
      rotating frame.

    Args:
      mu: Adimensional mass parameter of the system.
      
    Returns:
      L1, L2, L3, L4, L5: XY coordinates fo the Lagrange Points in the rotating 
        frame.
    """

    # L1 calculation
    coefL1 = [
        1, 4 * mu - 2, 6 * mu**2 - 6 * mu + 1,
        4 * mu**3 - 6 * mu**2 + 4 * mu - 1,
        mu**4 - 2 * mu**3 + 5 * mu**2 - 4 * mu + 2,
        2 * mu**3 - 3 * mu**2 + 3 * mu - 1
    ]
    L1roots = np.roots(coefL1)
    L1 = [0, 0]
    for i in L1roots:
        if np.isreal(i):
            L1 = np.array([i.real, 0])

    # L2 calculation
    coefL2 = [
        1, 4 * mu - 2, 6 * mu**2 - 6 * mu + 1,
        4 * mu**3 - 6 * mu**2 + 2 * mu - 1,
        mu**4 - 2 * mu**3 + mu**2 - 4 * mu + 2, 3 * mu - 3 * mu**2 - 1
    ]
    L2roots = np.roots(coefL2)
    L2 = [0, 0]
    for i in L2roots:
        if np.isreal(i):
            L2 = np.array([i.real, 0])

    # L3 calculation
    coefL3 = [
        1, 4 * mu - 2, 6 * mu**2 - 6 * mu + 1,
        4 * mu**3 - 6 * mu**2 + 2 * mu + 1,
        mu**4 - 2 * mu**3 + mu**2 + 4 * mu - 2, 3 * mu**2 - 3 * mu + 1
    ]
    L3roots = np.roots(coefL3)
    L3 = [0, 0]
    for i in L3roots:
        if np.isreal(i):
            L3 = np.array([i.real, 0])

    # L4 calculation
    L4 = np.array([1 / 2 - mu, 1 / 2 * m.sqrt(3)])

    # L5 calculation
    L5 = np.array([1 / 2 - mu, -1 / 2 * m.sqrt(3)])

    return np.array([L1, L2, L3, L4, L5])


def calc_jacobi(x, y, z, vx, vy, vz, mu):
    """Calculates Jacobi's Constant given the position and velocities.

    Args:
      x,y,z: Position coordinates of the Spacecraft.
      vx,vy,vz: Velocity of the spacecraft.
      mu: Adimensional mass parameter of the system.
      
    Returns:
      C: Jacobi's Constant.
    """

    r1 = distance('main', [x, y, z], mu)
    r2 = distance('secondary', [x, y, z], mu)
    C = x**2 + y**2 + 2 * (1 - mu) / r1 + 2 * mu / r2 - (vx**2 + vy**2 + vz**2)

    return C


def plot_3D_potential(mu, mainPos, secondaryPos, **kwargs):
    """Plots the 3D Potential of the 3BP in the rotational axis.

    Plots the 3D Potential (and optionally the Lagrange Points) of the 3BP in 
    the rotational axis

    Args:
      mu: Adimensional mass parameter of the system.
      mainPos: Position of the main, XYZ.
      secondaryPos: Position of the secondary, XYZ.
      **makeFigure (bool): Defaults to False. Boolean to make a separate figure 
        (True), or include it in the calling function's figure (False).
      **boolLagrange (bool): Defaults to True. Boolean to plot the Lagrange Points.
    """

    makeFigure = kwargs.get('makeFigure', False)
    boolLagrange = kwargs.get('boolLagrange', True)

    # prepare mesh for contourplot
    X, Y, Z = prepare_mesh_plot(bool3D=True)

    massiveObjectsHeight = np.nanmax(Z) + 0.5
    mainPos[2] = massiveObjectsHeight
    secondaryPos[2] = massiveObjectsHeight

    if makeFigure == True:
        fig, ax = prepare_plot_CR3BP(mainPos, secondaryPos, bool3D=True)

    if boolLagrange == True:
        LP = plot_lagrange(mu,
                           mainPos,
                           secondaryPos,
                           objectsHeight3D=massiveObjectsHeight,
                           axes=ax)

    surf = ax.plot_surface(X,
                           Y,
                           Z,
                           cmap=cm.gnuplot,
                           linewidth=0,
                           antialiased=False)  # contour for CjL
    ticks = np.linspace(np.nanmin(Z), np.nanmax(Z), 10)
    cbar = plt.colorbar(surf)
    cbar.set_label('- Jacobi Constant of ZVC')
    cbar.set_ticks(ticks)

    if makeFigure == True:
        ax.legend(loc='upper left')
        plt.title('Potential in the CR3BP')


def plot_zero_velocity_curves(mu, mainPos, secondaryPos, **kwargs):
    """Plots the ZVC of the 3BP in the rotational axis.

    Plots the ZVC (and optionally the Lagrange Points) of the 3BP in the 
    rotational axis.

    Args:
      mu: Adimensional mass parameter of the system.
      mainPos: Position of the main, XYZ.
      secondaryPos: Position of the secondary, XYZ.
      **makeFigure (bool): Defaults to False. Boolean to make a separate figure 
        (True), or include it in the calling function's figure (False).
      **boolLagrange (bool): Defaults to True. Boolean to plot the Lagrange Points.
      **scStateVector: Defaults to None. Initial coordinates of the Spacecraft 
        to calculate Jacobi's Constant.
    """

    makeFigure = kwargs.get('makeFigure', False)
    boolLagrange = kwargs.get('boolLagrange', True)
    scStateVector = kwargs.get('scStateVector', None)

    # prepare mesh for contourplot
    X, Y, Z = prepare_mesh_plot()

    if makeFigure == True:
        fig, ax = prepare_plot_CR3BP(mainPos, secondaryPos)

    if boolLagrange == True:
        LP = plot_lagrange(mu, mainPos, secondaryPos)

    if scStateVector == None:
        jacobiConst1 = calc_jacobi(LP[0, 0], LP[0, 1], 0, 0, 0, 0, mu)
        jacobiConst2 = calc_jacobi(LP[1, 0], LP[1, 1], 0, 0, 0, 0, mu)
        jacobiConst3 = calc_jacobi(LP[2, 0], LP[2, 1], 0, 0, 0, 0, mu)
        jacobiConst4 = calc_jacobi(LP[3, 0], LP[3, 1], 0, 0, 0, 0, mu)
        jacobiConst5 = calc_jacobi(LP[4, 0], LP[4, 1], 0, 0, 0, 0, mu)

        jacobiConstList = [jacobiConst3, jacobiConst2, jacobiConst1]

        plt.contour(X, Y, Z, jacobiConstList,
                    cmap=cm.gnuplot)  # contour for CjL

    else:
        jacobiConstSC = calc_jacobi(scStateVector[0], scStateVector[1],
                                    scStateVector[2], scStateVector[3],
                                    scStateVector[4], scStateVector[5], mu)
        Z[(jacobiConstSC <= Z[:, :]) | (Z[:, :] <
                                        (jacobiConstSC * 0.9))] = np.NaN

        levels = np.linspace(np.nanmin(Z), np.nanmax(Z), 100)
        ticks = np.linspace(np.nanmin(Z), np.nanmax(Z), 10)
        plt.contourf(X, Y, Z, cmap=cm.gnuplot, levels=levels)
        cbar = plt.colorbar()
        cbar.set_label('Jacobi Constant of ZVC')
        cbar.set_ticks(ticks)

    if makeFigure == True:
        ax.legend(loc='upper left')
        plt.title('Zero Velocity Curves')


def plot_lagrange(mu, mainPos, secondaryPos, **kwargs):
    """Plots the Lagrange Points in the 3BP in the rotational axis.

    Args:
      mu: Adimensional mass parameter of the system.
      mainPos: Position of the main, XYZ.
      secondaryPos: Position of the secondary, XYZ.
      **makeFigure (bool): Defaults to False. Boolean to make a separate figure 
        (True), or include it in the calling function's figure (False).
      **objectsHeight3D: Defaults to None. Z coordinate for the objects in 3D plot.
      **axes: Defaults to None. If 3D plot, the axes object needs to be passed.

    Returns:
      LP: Lagrange Points coordinates.
    """

    makeFigure = kwargs.get('makeFigure', False)
    objectsHeight3D = kwargs.get('objectsHeight3D', None)
    ax = kwargs.get('axes', None)

    LP = calculate_lagrange_points(mu)

    if makeFigure == True:
        fig, ax = prepare_plot_CR3BP(mainPos, secondaryPos)

    if objectsHeight3D == None:
        plt.scatter(LP[0, 0],
                    LP[0, 1],
                    marker='x',
                    color='r',
                    alpha=0.8,
                    zorder=9,
                    label='Lagrange Points')
        plt.scatter(LP[1, 0],
                    LP[1, 1],
                    marker='x',
                    color='r',
                    alpha=0.8,
                    zorder=9)
        plt.scatter(LP[2, 0],
                    LP[2, 1],
                    marker='x',
                    color='r',
                    alpha=0.8,
                    zorder=9)
        plt.scatter(LP[3, 0],
                    LP[3, 1],
                    marker='x',
                    color='r',
                    alpha=0.8,
                    zorder=9)
        plt.scatter(LP[4, 0],
                    LP[4, 1],
                    marker='x',
                    color='r',
                    alpha=0.8,
                    zorder=9)

    else:
        ax.scatter(LP[0, 0],
                   LP[0, 1],
                   objectsHeight3D,
                   s=50,
                   marker='x',
                   color='r',
                   alpha=0.8,
                   label='Lagrange Points')
        ax.scatter(LP[1, 0],
                   LP[1, 1],
                   objectsHeight3D,
                   s=50,
                   marker='x',
                   color='r',
                   alpha=0.8)
        ax.scatter(LP[2, 0],
                   LP[2, 1],
                   objectsHeight3D,
                   s=50,
                   marker='x',
                   color='r',
                   alpha=0.8)
        ax.scatter(LP[3, 0],
                   LP[3, 1],
                   objectsHeight3D,
                   s=50,
                   marker='x',
                   color='r',
                   alpha=0.8)
        ax.scatter(LP[4, 0],
                   LP[4, 1],
                   objectsHeight3D,
                   s=50,
                   marker='x',
                   color='r',
                   alpha=0.8)
        ax.text(LP[0, 0] - 0.05, LP[0, 1] - 0.15, objectsHeight3D, r'L$_1$')
        ax.text(LP[1, 0] + 0.05, LP[1, 1] - 0.05, objectsHeight3D, r'L$_2$')
        ax.text(LP[2, 0] - 0.2, LP[2, 1] - 0.05, objectsHeight3D, r'L$_3$')
        ax.text(LP[3, 0] - 0.05, LP[3, 1] + 0.05, objectsHeight3D, r'L$_4$')
        ax.text(LP[4, 0] - 0.05, LP[4, 1] - 0.15, objectsHeight3D, r'L$_5$')

    if makeFigure == True:
        ax.legend(loc='upper left')
        plt.title('Lagrange Points Position')

    return LP


def plot_CR3BP(mu, mainPos, secondaryPos, scStateVector, t_points, **kwargs):
    """Plots the CR3BP.

    Args:
      mu: Mass parameter of the system.
      referenceSystem: Reference System in which to plot, 'Rotating' or 
        'Inertial'. 
      mainPos: Position of the main, XYZ.
      secondarypos: Position of the secondary, XYZ.
      scStateVector: State vector of the SC, X Y Z VX VY VZ.
      t_points: Vector with the time points where the state vectors are.
      **boolZVC (bool): Defaults to False. Boolean to plot the ZVC and Lagrange 
        Points.
      **boolLagrange (bool): Defaults to False. Boolean to plot only the 
        Lagrange Points.
      **boolAnimation (list): List with information of the data as animation:

        - (bool): Defaults to False. True to create animation, otherwise False.

        - (bool): Defaults to False. True to save the animation, otherwise False.

        - Type of animation (str): Defaults to 'Standard'. Standard' or 
          'Comet'-like.

        - Reference System (str): Defaults to 'Rotating' or 'Inertial'.
    
    Returns:
        resultPlot: Plot object of the 3BP.
    """

    boolZVC = kwargs.get('boolZVC', False)
    boolLagrange = kwargs.get('boolLagrange', False)
    boolAnimation = kwargs.get('boolAnimation',
                               [False, False, 'Standard', 'Rotating'])

    fig, ax = prepare_plot_CR3BP(mainPos, secondaryPos)

    if boolZVC == True:
        plot_zero_velocity_curves(mu,
                                  mainPos,
                                  secondaryPos,
                                  scStateVector=[
                                      scStateVector[0, 0], scStateVector[1, 0],
                                      scStateVector[2, 0], scStateVector[3, 0],
                                      scStateVector[4, 0], scStateVector[5, 0]
                                  ])

    if boolLagrange == True:
        LP = plot_lagrange(mu, mainPos, secondaryPos)

    if boolAnimation[0] == True:
        mainPos = np.tile(mainPos, (t_points.shape[0], 1)).T
        secondaryPos = np.tile(secondaryPos, (t_points.shape[0], 1)).T

        if boolAnimation[3] == 'Inertial':
            mainPos = rotating_to_inertial_system(mainPos, adimAngVel,
                                                  t_points)
            secondaryPos = rotating_to_inertial_system(secondaryPos,
                                                       adimAngVel, t_points)
            scStateVector = rotating_to_inertial_system(
                scStateVector[:3], adimAngVel, t_points)

        resultPlot = plot_animation(fig,
                                    ax,
                                    mainPos,
                                    secondaryPos,
                                    scStateVector,
                                    t_points,
                                    animationType=boolAnimation[2],
                                    animationSave=boolAnimation[1])

    else:
        ax.plot(scStateVector[0], scStateVector[1], label='SC Trajectory')
        ax.scatter(mainPos[0],
                   mainPos[1],
                   marker='o',
                   s=40,
                   color='r',
                   alpha=0.7,
                   label='Primary')
        ax.scatter(secondaryPos[0],
                   secondaryPos[1],
                   marker='o',
                   s=20,
                   color='g',
                   alpha=0.7,
                   label='Secondary')
        finalTime, finalUnit = time_converter(t_points[-1] * adimTime,
                                              'seconds', 'days')
        plt.text(0.05,
                 0.05,
                 f't = {finalTime:.2f} {finalUnit}',
                 transform=ax.transAxes)
        ax.legend(loc=0, fontsize='x-small')
        plt.title('Orbit propagation for CR3BP (Earth - Moon)')
        resultPlot = ax

    return resultPlot


def prepare_plot_CR3BP(mainPos, secondaryPos, **kwargs):
    """Prepares the axis and fig part of the CR3BP plots.

    Args:
      mainPos: Position of the main, XYZ.
      secondarypos: Position of the secondary, XYZ.
      **bool3D (bool): Defaults to False. Boolean to prepare for the 3D plot 
        instead of the planar.
    
    Returns:
      fig: Fig object of the plot.
      ax: Axes object of the plot.
    """

    bool3D = kwargs.get('bool3D', False)

    if bool3D == True:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax = plt.gca(projection='3d')

    else:
        fig, ax = plt.subplots()
        ax = plt.gca()

    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_aspect('equal')

    return fig, ax


def prepare_mesh_plot(**kwargs):
    """Prepares the mesh for the potential calculations in the CR3BP plots.
      
    Args:
      **bool3D (bool): Defaults to False. Boolean to prepare for the 3D plot 
        instead of the planar.
  
    Returns:
      X, Y, Z: Values for the mesh of the potential.
    """

    bool3D = kwargs.get('bool3D', False)

    # prepare mesh for contourplot
    delta = 0.001
    x = np.arange(-1.5, 1.5, delta)
    y = np.arange(-1.5, 1.5, delta)
    X, Y = np.meshgrid(x, y)
    if bool3D == False:
        Z = calc_jacobi(X, Y, 0, 0, 0, 0, mu)
    else:
        Z = -calc_jacobi(X, Y, 0, 0, 0, 0, mu)
        Z[Z[:, :] < (-7)] = -7

    return X, Y, Z


def plot_animation(fig, ax, mainPos, secondaryPos, scStateVector, t_points,
                   **kwargs):
    """Plots the animation of the 3BP and optionally saves it.
    
    Args:
      fig: Fig object of the plot.
      ax: Axes object of the plot.
      mainPos: Position of the main, XYZ.
      secondarypos: Position of the secondary, XYZ.
      scStateVector: State vector of the SC, X Y Z VX VY VZ.
      t_points: Vector with the time points where the state vectors are.
      **animationType (str); Defaults to 'Standard'. 'Standard' or 'Comet'-like.
      **animationSave (bool): Defaults to False. Boolean to save or not the 
        animation.

    Returns:
      resultPlot: Animated plot object of the 3BP.
    """

    animationType = kwargs.get('animationType', 'Standard')
    animationSave = kwargs.get('animationSave', 'False')

    mainPosition, = ax.plot([], [],
                            marker='o',
                            markersize=10,
                            color='r',
                            alpha=1,
                            label='Primary')
    secondaryPosition, = ax.plot([], [],
                                 marker='o',
                                 markersize=5,
                                 color='g',
                                 alpha=1,
                                 label='Secondary')
    scPositionHistory, = ax.plot([], [], c='b', zorder=9)
    scPositionCurrent, = ax.plot([], [],
                                 c='b',
                                 marker='o',
                                 markersize=2,
                                 markerfacecolor='tab:orange',
                                 markeredgecolor='tab:orange',
                                 label='SC Trajectory',
                                 zorder=10)
    time_points = ax.text(0.05, 0.05, '', transform=ax.transAxes)
    ax.legend(loc='upper left')
    plt.title('Orbit propagation for 3BP (Earth - Moon)')
    timePoints, unit = time_converter(t_points * adimTime, 'seconds', 'days')

    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=100,
                    metadata=dict(artist='Me'),
                    bitrate=1800,
                    extra_args=['-vcodec', 'libx264'])

    def init():
        """initialize animation"""
        scPositionCurrent.set_data([], [])
        scPositionHistory.set_data([], [])
        mainPosition.set_data([], [])
        secondaryPosition.set_data([], [])
        time_points.set_text('')

        return scPositionCurrent, scPositionHistory, mainPosition, secondaryPosition, time_points

    # animation function.  This will be called sequentially with the frame number
    def animate(i):
        time_points.set_text(f't = {timePoints[i]:.2f} {unit}')
        mainPosition.set_data(mainPos[0, i], mainPos[1, i])
        secondaryPosition.set_data(secondaryPos[0, i], secondaryPos[1, i])
        if animationType == 'Standard':
            scPositionCurrent.set_data(scStateVector[0, i],
                                       scStateVector[1, i])
            scPositionHistory.set_data(scStateVector[0, :i],
                                       scStateVector[1, :i])
        elif animationType == 'Comet':
            scPositionCurrent.set_data(scStateVector[0, i],
                                       scStateVector[1, i])
            if i < 500:
                scPositionHistory.set_data(scStateVector[0, :i],
                                           scStateVector[1, :i])
            else:
                scPositionHistory.set_data(scStateVector[0, i - 500:i],
                                           scStateVector[1, i - 500:i])
        return scPositionCurrent, scPositionHistory, mainPosition, secondaryPosition, time_points

    resultPlot = animation.FuncAnimation(fig,
                                         animate,
                                         init_func=init,
                                         frames=len(scStateVector[0]),
                                         interval=10,
                                         blit=True)

    if animationSave == True:
        # Save as mp4. This requires mplayer or ffmpeg to be installed
        resultPlot.save('3bp.mp4', writer=writer)

    return resultPlot


def inertial_to_rotating_system(objectPositionInertial, angVel_system, t):
    """Transforms the coordinates of an object in the CR3BP from the inertial 
    frame to the rotating.
  
    Args:
      objectPositionInertial: Coordinates of the object in the inertial frame, 
        XYZ and time.
      angVel_system: Angular velocity of the 2-body system.
      t: Time where the coordinates are.
  
    Returns:
      objectPositionRotating: Coordinates of the object in the rotating frame, 
        XYZ and time.
    """

    objectPositionRotating = np.empty(
        (objectPositionInertial.shape[0], t.shape[0]))

    for timePoint in range(0, t.shape[0]):
        C = np.array([[
            np.cos(angVel_system * t[timePoint]),
            np.sin(angVel_system * t[timePoint]), 0
        ],
                      [
                          -np.sin(angVel_system * t[timePoint]),
                          np.cos(angVel_system * t[timePoint]), 0
                      ], [0, 0, 1]])
        objectPositionRotating[:3, timePoint] = np.dot(
            C, objectPositionInertial[:3, timePoint])[:]
        if objectPositionInertial.shape[0] == 6:
            objectPositionRotating[3:6, timePoint] = np.dot(
                np.dot(
                    np.array([[0, angVel_system, 0], [-angVel_system, 0, 0],
                              [0, 0, 0]]), C),
                objectPositionInertial[:3, timePoint])[:] + np.dot(
                    C, objectPositionInertial[3:6, timePoint])[:]

    return objectPositionRotating


def rotating_to_inertial_system(objectPositionRotating, angVel_system, t):
    """Transforms the coordinates of an object in the CR3BP from the rotating 
    frame to the inertial.
  
    Args:
      objectPositionRotating: Coordinates of the object in the rotating frame, 
        XYZ and time.
      angVel_system: Angular velocity of the 2-body system.
      t: Time where the coordinates are.
  
    Returns:
      objectPositionInertial: Coordinates of the object in the inertial frame, 
        XYZ and time.
    """

    objectPositionInertial = np.empty(
        (objectPositionRotating.shape[0], t.shape[0]))

    for timePoint in range(0, t.shape[0]):
        C = np.array([[
            np.cos(angVel_system * t[timePoint]),
            np.sin(angVel_system * t[timePoint]), 0
        ],
                      [
                          -np.sin(angVel_system * t[timePoint]),
                          np.cos(angVel_system * t[timePoint]), 0
                      ], [0, 0, 1]])
        objectPositionInertial[:3, timePoint] = np.dot(
            C.T, objectPositionRotating[:3, timePoint])[:]
        if objectPositionRotating.shape[0] == 6:
            objectPositionInertial[3:6, timePoint] = np.dot(
                np.dot(
                    np.array([[0, -angVel_system, 0], [angVel_system, 0, 0],
                              [0, 0, 0]]), C.T),
                objectPositionRotating[:3, timePoint])[:] + np.dot(
                    C.T, objectPositionRotating[3:6, timePoint])[:]

    return objectPositionInertial


def time_converter(time, currentUnit, objectiveUnit):
    """Converts time from one unit to another.

    Converts time values from the time variable that are in currentUnit to the 
    unit in objectiveUnit.
  
    Args:
      time: Time value or list of time values.
      currentUnit (str): Unit in which the values of time are expressed 
        (seconds, hours, days, years).
      objectiveUnit (str): Unit to which to convert the time values (seconds, 
        hours, days, years).
  
    Returns:
      resultTime: Time value or list of time values in the desired unit.
      objectiveUnit (str): Unit in which the times in resultTime are expressed 
        (seconds, hours, days, years).
    """

    hour_in_seconds = 3600
    day_in_hours = 24
    year_in_days = 365

    if currentUnit == 'seconds':
        if objectiveUnit == 'hours':
            resultTime = time / hour_in_seconds

        elif objectiveUnit == 'days':
            resultTime = time / (hour_in_seconds * day_in_hours)

        elif objectiveUnit == 'years':
            resultTime = time / (hour_in_seconds * day_in_hours * year_in_days)

    elif currentUnit == 'hours':
        if objectiveUnit == 'seconds':
            resultTime = time * hour_in_seconds

        elif objectiveUnit == 'days':
            resultTime = time / (day_in_hours)

        elif objectiveUnit == 'years':
            resultTime = time / (day_in_hours * year_in_days)

    elif currentUnit == 'days':
        if objectiveUnit == 'seconds':
            resultTime = time * day_in_hours * hour_in_seconds

        elif objectiveUnit == 'hours':
            resultTime = time * day_in_hours

        elif objectiveUnit == 'years':
            resultTime = time / year_in_days

    elif currentUnit == 'years':
        if objectiveUnit == 'seconds':
            resultTime = time * year_in_days * day_in_hours * hour_in_seconds
        elif objectiveUnit == 'hours':
            resultTime = time * year_in_days * day_in_hours
        elif objectiveUnit == 'days':
            resultTime = time * year_in_days

    return [resultTime, objectiveUnit]
