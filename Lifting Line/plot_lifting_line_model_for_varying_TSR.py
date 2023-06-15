# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 16:07:52 2023

@author: MFekry

This script creates all plots for TSR=6, 8, 10
"""
import numpy as np
import matplotlib.pyplot as plt
from lifting_line_model_functions import create_rotor_geometry, solve_lifting_line_system_matrix_approach,changewakeconvectionspeedplot

plt.close('all')

# Generate normal wake geometry (control wake length and the element number and convection of the wake)
NELEMENTS = 10
Nrotations = 25
Rotation_resoltion = 20
aw = 0.0

s_Array = np.linspace(0, np.pi, NELEMENTS)
# s_len_Array = np.linspace(0,1, NELEMENTS)
print(s_Array)

s_cos_array = []
for i in range(len(s_Array)):
    s_cos_array.append((-0.8/2.0 * (np.cos(s_Array[i])-1)+0.2) * 50.0 )

print(s_cos_array)
# plt.plot(s_len_Array,s_cos_array)
maxradius = max(s_cos_array)
theta_Array = np.linspace(0, Nrotations*2*np.pi, Rotation_resoltion)

print(theta_Array)

def return_lift_solution(TSR):
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

# Plot 'Ftan vs r_R'
lift_solution_6, nonDim_Gamma_6, nonDim_Force_6 = return_lift_solution(TSR=6)
lift_solution_8, nonDim_Gamma_8, nonDim_Force_8 = return_lift_solution(TSR=8)
# lift_solution_10, nonDim_Gamma_10, nonDim_Force_10 = return_lift_solution(TSR=10)

x = lift_solution_6['r_R']
y_6 = lift_solution_6['Ftan']/nonDim_Force_6
y_8 = lift_solution_8['Ftan']/nonDim_Force_8
# y_10 = lift_solution_10['Ftan']/nonDim_Force_10

plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')

plt.ylabel('Ftan',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('Ftan vs r_R.png',dpi=600,bbox_inches ='tight')

# Plot 'Fnorm vs r_R'
y_6 = lift_solution_6['Fnorm']/nonDim_Force_6
y_8 = lift_solution_8['Fnorm']/nonDim_Force_8
# y_10 = lift_solution_10['Fnorm']/nonDim_Force_10

plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')
plt.ylabel('Fnorm',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('Fnorm vs r_R.png',dpi=600,bbox_inches ='tight')
# plt.close()

# Plot 'Gamma vs r_R'

y_6 = lift_solution_6['Gamma']/nonDim_Gamma_6
y_8 = lift_solution_8['Gamma']/nonDim_Gamma_8
# y_10 = lift_solution_10['Gamma']/nonDim_Gamma_10

plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')
plt.ylabel('Gamma',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('Gamma vs r_R.png',dpi=600,bbox_inches ='tight')

# Plot 'a vs r_R'
y_6 = lift_solution_6['a']
y_8 = lift_solution_8['a']
# y_10 = lift_solution_10['a']
plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')

plt.ylabel('a',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('a vs r_R.png',dpi=600,bbox_inches ='tight')
# plt.close()

# 'aline vs r_R'
y_6 = lift_solution_6['aline']
y_8 = lift_solution_8['aline']
# y_10 = lift_solution_10['aline']
plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')

plt.ylabel('aline',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('aline vs r_R.png',dpi=600,bbox_inches ='tight')
# plt.close()

# 'CP_array vs r_R'
y_6 = lift_solution_6['CP_array']
y_8 = lift_solution_8['CP_array']
# y_10 = lift_solution_10['CP_array']
plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')
plt.ylabel('CP',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('CP_array vs r_R.png',dpi=600,bbox_inches ='tight')
# plt.close()

# 'CT_array vs r_R'
y_6 = lift_solution_6['CT_array']
y_8 = lift_solution_8['CT_array']
# y_10 = lift_solution_10['CT_array']
plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')
plt.ylabel('CT',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('CT_array vs r_R.png',dpi=600,bbox_inches ='tight')
# plt.close()

# 'Alpha vs r_R'
y_6 = lift_solution_6['alpha']
y_8 = lift_solution_8['alpha']
# y_10 = lift_solution_10['alpha']
plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')
plt.ylabel('Alpha (deg)',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('Alpha vs r_R.png',dpi=600,bbox_inches ='tight')
# plt.close()

# 'inflow vs r_R'
y_6 = np.array(lift_solution_6['inflow']) * 180 / np.pi
y_8 = np.array(lift_solution_8['inflow']) * 180 / np.pi
# y_10 = lift_solution_10['inflow']
plt.figure(figsize=(width,length))
plt.grid(axis='both')
plt.plot(x, y_6, linewidth=2, linestyle="-", alpha=1, label='TSR=6')
plt.plot(x, y_8, linewidth=2, linestyle="-", alpha=1, label='TSR=8')
# plt.plot(x, y_10, linewidth=2, linestyle="-", alpha=1, label='TSR=10')
plt.ylabel('inflow angle (deg)',fontsize=labelFontsize)
plt.xlabel('r_R',fontsize=labelFontsize)
plt.tick_params(axis='both', which='both', labelsize=axlabelsize)
plt.legend(fontsize=labelFontsize)
plt.savefig('inflow vs r_R.png',dpi=600,bbox_inches ='tight')
# plt.close()