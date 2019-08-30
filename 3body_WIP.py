# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 10:34:29 2019

@author: Roger Gutierrez Ramon
"""

# Import general packages and modules
import numpy as np
import math as m

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
from matplotlib import animation
from matplotlib.animation import FFMpegWriter
from scipy.integrate import solve_ivp
from scipy import optimize

from astropy import constants as const
import spiceypy as spice

# Import personal packages and modules
from astro_equations.variables_3bp import *
from astro_equations.support_functions_3bp import *

# Set general animation\figure parameters
matplotlib.rcParams['figure.dpi'] = 200
matplotlib.rcParams['animation.bitrate'] = 1800
matplotlib.rcParams['figure.figsize'] = 10, 10

# Define initial conditions for the ode
Y0, timeSpan, timePoints, boolZVC, boolLagrange, boolAnimation = define_initial_conditions(
    mu)

# Integration
orbit = solve_ivp(state_vector_CR3BP,
                  timeSpan,
                  Y0,
                  method='RK45',
                  t_eval=timePoints,
                  rtol=1e-12,
                  atol=1e-12)

# Unpacking of variables from ode and adapting primary and secondary objects to
# time sequence
rot_SC = np.array(
    [orbit.y[0], orbit.y[1], orbit.y[2], orbit.y[3], orbit.y[4], orbit.y[5]])
t_points = np.array(orbit.t)

# Plotting results and wrapping up
#plot3DPotential(mu, rot_PrimaryPos, rot_SecondaryPos, makeFigure=True)
#plotZVC(mu, rot_PrimaryPos, rot_SecondaryPos, makeFigure=True)
#resultPlot = plotCR3BProt(mu, rot_PrimaryPos, rot_SecondaryPos, rot_SC, t_points, boolZVC=True, boolAnimation=[True,'Standard'])
resultPlot = plot_CR3BP(mu,
                       rot_PrimaryPos,
                       rot_SecondaryPos,
                       rot_SC,
                       t_points,
                       boolLagrange=boolLagrange,
                       boolZVC=boolZVC,
                       boolAnimation=boolAnimation)

plt.show()