# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 16:07:52 2023

@author: MFekry

This script creates all plots for TSR=6, 8, 10
"""
import numpy as np
import matplotlib.pyplot as plt
from lifting_line_model_functions import create_rotor_geometry, solve_lifting_line_system_matrix_approach,changewakeconvectionspeedplot


def return_lift_solution(NELEMENTS, Nrotations, Rotation_resolution):
    # Generate normal wake geometry (control wake length and the element number and convection of the wake)
    aw = 0.0
    TSR = 6

    s_Array = np.linspace(0, np.pi, NELEMENTS)
    # s_len_Array = np.linspace(0,1, NELEMENTS)
    print(s_Array)

    s_cos_array = []
    for i in range(len(s_Array)):
        s_cos_array.append((-0.8/2.0 * (np.cos(s_Array[i])-1)+0.2) * 50.0 )

    print(s_cos_array)
    # plt.plot(s_len_Array,s_cos_array)
    maxradius = max(s_cos_array)
    theta_Array = np.linspace(0, Nrotations * 2 * np.pi, Rotation_resolution)

    print(theta_Array)

    Uinf = 1
    nblades = 1

    rotor_wake_system = create_rotor_geometry(s_cos_array, maxradius, TSR/(1-aw), Uinf, theta_Array,nblades)

    wind = [Uinf,0 , 0]
    Omega = TSR* Uinf/maxradius
    AirDen = 1.225
    lift_solution = solve_lifting_line_system_matrix_approach(rotor_wake_system, wind, Omega,maxradius,nblades)

    nonDim_Force = 0.5* AirDen* Uinf**2 *maxradius
    nonDim_Gamma = Uinf**2*np.pi/(nblades*Omega)

    return lift_solution, nonDim_Gamma, nonDim_Force

scale=1
width=9*scale
length=6*scale
labelFontsize=16
axlabelsize=16

NELEMENTS = 10
Nrotations = 25
Rotation_resolution = 20

# Plot 'a vs r_R' for varying number of elements
plt.figure(figsize=(width,length))
plt.grid(axis='both')
elements = [6, 8, 10]

for n in elements:
    lift_solution, nonDimGamma, nonDimForce = return_lift_solution(
        NELEMENTS=n,
        Nrotations=Nrotations,
        Rotation_resolution=Rotation_resolution
    )

    plt.plot(
        lift_solution['r_R'],
        lift_solution['a'],
        linestyle='-',
        alpha=1,
        label=f'{n} elements'
    )

plt.ylabel('a',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('SENS_ANALYSIS_a vs r_R varying elements.png', dpi=600, bbox_inches='tight')
# plt.close()


# Plot 'a vs r_R' for varying number of rotations
plt.figure(figsize=(width,length))
plt.grid(axis='both')
rotations = [15, 25, 40]

for n in rotations:
    lift_solution, nonDimGamma, nonDimForce = return_lift_solution(
        NELEMENTS=NELEMENTS,
        Nrotations=n,
        Rotation_resolution=Rotation_resolution
    )

    plt.plot(
        lift_solution['r_R'],
        lift_solution['a'],
        linestyle='-',
        alpha=1,
        label=f'{n} rotations'
    )

plt.ylabel('a',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('SENS_ANALYSIS_a vs r_R varying rotations.png', dpi=600, bbox_inches='tight')
# plt.close()

