# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:08:28 2023

@author: MFekry
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os



def create_airfoil_geometry(a_array, Flat_blade_chord, NELEMENTS, alpha, DXW, SX, SZ):
  SN = np.sin(alpha)
  CS = np.cos(alpha)  
  filaments = []
  ring = []
  controlpoints = []
  vortexpoints = []
  shedwakepoints =[]

  # for krot in range(nblades): 

  #   angle_rotation = 2*np.pi/nblades*krot
  #   cosrot = np.cos(angle_rotation)
  #   sinrot = np.sin(angle_rotation)

  for i in range(len(a_array)-1):
      # print('i = ',i)
      c = (a_array[i] + a_array[i+1])/2;
      delta_c = (a_array[i+1] - a_array[i])
      c_center = c/delta_c
      # define controlpoints     
      temp1 = {'coordinates': [ (c_center + 0.25)*CS + SX ,  0  , -1*(c_center + 0.25)*SN + SZ ] , 'ChortCenter': c_center , 'DeltaChord': delta_c  }
      controlpoints.append( temp1 );
      # define vortexpoints
      temp1 = {'coordinates': [(c_center - 0.25)*CS + SX ,  0  , -1*(c_center - 0.25)*SN + SZ ] , 'strenght': 0.0 }
      vortexpoints.append( temp1 );
      # # define shedwakepoints
      # temp1 = {'coordinates': [(c_center + DXW)*CS + SX ,  0  , -1*(c_center + DXW)*SN + SZ ] }
      # shedwakepoints.append( temp1 );
  
  return ({'controlpoints': controlpoints ,  'vortexpoints':vortexpoints })#,  'shedwakepoints':shedwakepoints})
 
    
 
def VOR2D(Gamma,x,z,xi,zi):
    rz = z-zi
    rx = x-xi
    r = np.sqrt(rx**2 + rz**2)
    if r < 0.001:
        r = 0.001
    Cir = 0.5*Gamma/(np.pi * r)
    U = Cir * (rz/r)
    W = -Cir * (rx/r)
    return [U,W]


def DownWash_CP(NELEMENTS, x, z, Wake_Geometry, It1, It2):  
    U = 0
    W = 0
    # print('Wake_Geometry = ', Wake_Geometry)
    for It in range(It1,It2):
        # print('It - 1 = ', It)
        Ucp = 0
        Wcp = 0
        xi = Wake_Geometry['coordinates'][0]
        zi = Wake_Geometry['coordinates'][2]
        # wait = input("Press the <ENTER> key to continue...")
        for icp in range(NELEMENTS):
            Gamma1 = Wake_Geometry['strenght']
            U1,W1 = VOR2D(Gamma1,x[icp],z[icp],xi,zi)
            Ucp = Ucp + U1
            Wcp = Wcp + W1
        U = U + Ucp
        W = W + Wcp
    return [U,W]

def DownWash_WR(x, z, Wake_Geometry, It1, It2):  
    U = 0
    W = 0
    # print('Wake_Geometry = ', Wake_Geometry)
    for It in range(It1,It2):
        # print('It - 1 = ', It)
        Ucp = 0
        Wcp = 0
        # wait = input("Press the <ENTER> key to continue...")
        for IWG in range(len(Wake_Geometry)):
            xi = Wake_Geometry[IWG]['coordinates'][0]
            zi = Wake_Geometry[IWG]['coordinates'][2]
            Gamma1 = Wake_Geometry[IWG]['strenght']
            U1,W1 = VOR2D(Gamma1,x,z,xi,zi)
            Ucp = Ucp + U1
            Wcp = Wcp + W1
        U = U + Ucp
        W = W + Wcp
    return [U,W]
 
def time2semichord(time, Uinf, chord):
    return 2*Uinf*time/chord

def semichord2time(s, Uinf, chord):
    return s/2/Uinf*chord

# determining X and Y terms for recursive marching formula for approximation of Duhamel's integral 
def duhamel_approx(Xi,Yi,delta_s,delta_alpha,order=2,A1=0.3,A2=0.7,b1=0.14,b2=0.53):
    # A1=0.165,A2=0.335,b1=0.0455,b2=0.3
    # determine the next values of X and Y, named Xip1 and Yip1
    if order==1:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha
    elif order==2:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha*np.exp(-b1*delta_s/2)
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha*np.exp(-b2*delta_s/2)        
    else:
        Xip1= Xi*np.exp(-b1*delta_s)+A1*delta_alpha*((1+4*np.exp(-b1*delta_s/2)+np.exp(-b1*delta_s))/6)
        Yip1= Yi*np.exp(-b2*delta_s)+A2*delta_alpha*((1+4*np.exp(-b2*delta_s/2)+np.exp(-b2*delta_s))/6)
    
    return Xip1,Yip1


# define function for circulatory force, potential flow
def circulatory_normal_force(dCn_dalpha,alpha_equivalent,alpha0):
    return dCn_dalpha*(alpha_equivalent-alpha0)


def vorticity_to_velocity_field(x, z, gammas, vortex_x_coordinates, vortex_z_coordinates):

    Vx = sum(
            gammas * (z - vortex_z_coordinates) /
            (2 * np.pi * ((x - vortex_x_coordinates)**2 + (z - vortex_z_coordinates)**2))
    )

    Vz = sum(
            -1 * gammas * (x - vortex_x_coordinates) /
            (2 * np.pi * ((x - vortex_x_coordinates)**2 + (z - vortex_z_coordinates)**2))
    )

    return Vx, Vz

# if __name__ == "__main__": 
#     NELEMENTS =5
#     # NELEMENTS += 1
#     Flat_blade_chord = 60 
#     Flat_blade_Leading_edge = 0.0
#     airfoil_Array = np.linspace(0, 1*np.pi, NELEMENTS)
#     airfoil_cos_array = []
#     Uinf = 1 #m/s
#     alpha = 90*np.pi/180 # angle of attack
    
#     for i in range(len(airfoil_Array)+1):
#         # airfoil_cos_array.append((-1/2.0* (np.cos(airfoil_Array[i])-1) + Flat_blade_Leading_edge) * Flat_blade_chord )
#         airfoil_cos_array.append((i)*Flat_blade_chord/NELEMENTS + Flat_blade_Leading_edge)
#     # print(s_cos_array)
#     # plt.plot(s_Array,s_cos_array)
    
#     Airfoil_geometry = create_airfoil_geometry(airfoil_cos_array,Flat_blade_chord,NELEMENTS)
    
#     Matrix = np.zeros([len(Airfoil_geometry['controlpoints']),len(Airfoil_geometry['vortexpoints'])])
#     MatrixU = np.zeros([len(Airfoil_geometry['controlpoints']),len(Airfoil_geometry['vortexpoints'])])
#     MatrixW = np.zeros([len(Airfoil_geometry['controlpoints']),len(Airfoil_geometry['vortexpoints'])])
#     Gamma = np.ones([len(Airfoil_geometry['controlpoints']),1])
#     Gammanew = Gamma
#     # Matrix = np.zeros([5,5])
    
#     for icp in range(len(Airfoil_geometry['controlpoints'])):
#         delta_c = 1
#         for jvortex in range(len(Airfoil_geometry['vortexpoints'])):
#             # r = (Airfoil_geometry['controlpoints'][icp]['coordinates'][0] - Airfoil_geometry['vortexpoints'][jvortex]['coordinates'][0])**2
#             # Matrix[icp][jvortex] = -0.5 /(r) * (Airfoil_geometry['controlpoints'][icp]['coordinates'][0] - Airfoil_geometry['vortexpoints'][jvortex]['coordinates'][0])
#             UW = VOR2D(Gamma[jvortex],Airfoil_geometry['controlpoints'][icp]['coordinates'][0],0,Airfoil_geometry['vortexpoints'][jvortex]['coordinates'][0],0)
#             MatrixU[icp][jvortex] = UW[0]
#             MatrixW[icp][jvortex] = UW[1]
            
#             # 
#     print('a = ', MatrixW)
#     RHS = Gamma
#     GammaNew = -1*Uinf*delta_c*np.sin(alpha)*np.matmul(np.linalg.inv(MatrixW), RHS)
#     print('Gamma = ', GammaNew)
    
    
    














