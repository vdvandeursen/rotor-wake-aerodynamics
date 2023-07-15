# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 23:27:08 2023

@author: MFekry


Ref: Aerodynamics - Katz & Plotkin appendix D program 15 and 16 with (sec 13.10)
"""
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# import system
from Unsteady_Aerodynamics_function import create_airfoil_geometry,VOR2D, DownWash_CP, DownWash_WR, time2semichord, semichord2time, duhamel_approx, circulatory_normal_force, vorticity_to_velocity_field
plt.close('all')

if __name__ == "__main__": 
    NELEMENTS = 5
    # ShedWake = 2
    # NELEMENTS += 1
    Flat_blade_chord = 50 
    Flat_blade_Leading_edge = 0.0
    airfoil_Array = np.linspace(0, 1*np.pi, NELEMENTS)
    
    Airfoil_geometry =[]
    Nsteps = 80
    rho = 1.225
    
    # Variable declaration
    Gamma = np.ones([NELEMENTS+1,Nsteps])
    Gammanew = Gamma
    RHS = np.zeros([NELEMENTS+1,Nsteps])
    MatrixU = np.zeros([Nsteps,NELEMENTS+1,NELEMENTS+1])
    MatrixW = np.zeros([Nsteps,NELEMENTS+1,NELEMENTS+1])
    DownWashU = np.zeros([Nsteps,1])
    DownWashW = np.zeros([Nsteps,1])
    WWake = np.zeros([Nsteps,NELEMENTS])
    Qt_upper = np.zeros([Nsteps,NELEMENTS])
    Qt_lower = np.zeros([Nsteps,NELEMENTS])
    DGamma = np.zeros([Nsteps,NELEMENTS])
    delta_P = np.zeros([Nsteps,NELEMENTS])
    Drag_panel = np.zeros([Nsteps,NELEMENTS])
    Lift = np.zeros([Nsteps,1])
    Drag = np.zeros([Nsteps,1])
    CL = np.zeros([Nsteps,1])
    CD = np.zeros([Nsteps,1])
    Wake_Geometry = []
    
    # Wind Conditions
    Uinf = 50 #m/s
    alpha = 5*np.pi/180 # angle of attack
    SN = np.sin(alpha)
    CS = np.cos(alpha)
    
    # time step
    Dt = 0.25 * Flat_blade_chord/Uinf/NELEMENTS
    DXW = 0.3 * Uinf * Dt
        
    for It in range(Nsteps):
        T = It * Dt
        SX = -Uinf*T
        SZ = 0.0
        airfoil_chord_array = []
        xcp =np.zeros([NELEMENTS,1])
        zcp =np.zeros([NELEMENTS,1])
        xvx =np.zeros([NELEMENTS,1])
        zvx =np.zeros([NELEMENTS,1])
        
        # create chord geometery
        for i in range(len(airfoil_Array)+1):
            # airfoil_chord_array.append((-1/2.0* (np.cos(airfoil_Array[i])-1) + Flat_blade_Leading_edge) * Flat_blade_chord )
            airfoil_chord_array.append((i)*Flat_blade_chord/NELEMENTS + Flat_blade_Leading_edge)
        Airfoil_geometry.append(create_airfoil_geometry(airfoil_chord_array,Flat_blade_chord,NELEMENTS, alpha, DXW, SX, SZ))
        # Matrix = np.zeros([Nsteps, len(Airfoil_geometry['controlpoints']),len(Airfoil_geometry['vortexpoints'])])
        
        
        # shad wake properties
        delta_c = 1
        Wake_Vort_X =  (NELEMENTS*delta_c + DXW)* CS + SX
        Wake_Vort_Z = -(NELEMENTS*delta_c + DXW)* SN + SZ
        
        temp1 = {'coordinates': [Wake_Vort_X,  0  , Wake_Vort_Z ] , 'strenght': 1.0}
        Wake_Geometry.append(temp1)
       
        
        # Calculate matrix aij coefficients
        
        print('Step = ', It)
        for icp in range(NELEMENTS):
            xcp[icp] = Airfoil_geometry[It]['controlpoints'][icp]['coordinates'][0]
            zcp[icp] = Airfoil_geometry[It]['controlpoints'][icp]['coordinates'][2]
            
            for jvortex in range(NELEMENTS):
                # r = (Airfoil_geometry['controlpoints'][icp]['coordinates'][0] - Airfoil_geometry['vortexpoints'][jvortex]['coordinates'][0])**2
                # Matrix[icp][jvortex] = -0.5 /(r) * (Airfoil_geometry['controlpoints'][icp]['coordinates'][0] - Airfoil_geometry['vortexpoints'][jvortex]['coordinates'][0])
                UW = VOR2D(Gamma[jvortex][It],Airfoil_geometry[It]['controlpoints'][icp]['coordinates'][0],Airfoil_geometry[It]['controlpoints'][icp]['coordinates'][2],Airfoil_geometry[It]['vortexpoints'][jvortex]['coordinates'][0],Airfoil_geometry[It]['vortexpoints'][jvortex]['coordinates'][2])
                # print('UW = ' ,UW)
                MatrixU[It][icp][jvortex] = UW[0]
                MatrixW[It][icp][jvortex] = UW[1]
                MatrixU[It][NELEMENTS][jvortex] = 1
                MatrixW[It][NELEMENTS][jvortex] = 1
                
                
               
                # UW = DownWash(Airfoil_geometry[It]['controlpoints'][icp]['coordinates'][0],Airfoil_geometry[It]['controlpoints'][icp]['coordinates'][2],0,It-1,Gamma[jvortex][It],Airfoil_geometry[It]['shedwakepoints'][jvortex]['coordinates'][0],Airfoil_geometry[It]['shedwakepoints'][jvortex]['coordinates'][2]) 
            UW = VOR2D(Gamma[NELEMENTS][It],Airfoil_geometry[It]['controlpoints'][icp]['coordinates'][0],Airfoil_geometry[It]['controlpoints'][icp]['coordinates'][2],Wake_Geometry[It]['coordinates'][0],Wake_Geometry[It]['coordinates'][2])
            MatrixU[It][icp][NELEMENTS] = UW[0]
            MatrixW[It][icp][NELEMENTS] = UW[1]
            
            if It>0:
                UW = DownWash_CP(NELEMENTS, xcp,zcp, Wake_Geometry[It], 0, It)
                WWake[It,icp] = UW[0]*SN + UW[1]*CS  
             
                
            # RHS vector calcualtions 
            RHS[icp][It] = -Uinf*SN - WWake[It,icp]
            
        RHS[NELEMENTS][It] -=   Gamma[NELEMENTS,It]
        MatrixU[It][NELEMENTS][NELEMENTS] = 1
        MatrixW[It][NELEMENTS][NELEMENTS] = 1 
        Gamma[:,It] = np.matmul(np.linalg.inv(MatrixW[It][:][:]), RHS[:,It])
        Wake_Geometry[It]['strenght'] = Gamma[-1,It]
        
        for jvortex in range(NELEMENTS):
            Airfoil_geometry[It]['vortexpoints'][jvortex]['strenght'] =  Gamma[jvortex,It]
        
        
        # Wake Rollup
        if It > 0 :
            UWT = [0,0]
            UWvx = np.zeros([It,2])
            UWWR = np.zeros([It,2])
            for IWR in range(It):
                for jvortex in range(NELEMENTS):
                    xvx[jvortex] = Airfoil_geometry[IWR]['vortexpoints'][jvortex]['coordinates'][0]
                    zvx[jvortex] = Airfoil_geometry[IWR]['vortexpoints'][jvortex]['coordinates'][2]
                    
                    UW = VOR2D(Gamma[jvortex][IWR],Wake_Geometry[IWR]['coordinates'][0],Wake_Geometry[IWR]['coordinates'][2],xvx[jvortex],zvx[jvortex])
                    UW1 = DownWash_WR(Wake_Geometry[IWR]['coordinates'][0],Wake_Geometry[IWR]['coordinates'][2], Wake_Geometry, 0, It)
                    
                    UWvx[IWR,0] +=  (UW[0] + UW1[0]) # U
                    UWvx[IWR,1] +=  (UW[1] + UW1[1]) # W
                    
               
                Wake_Geometry[IWR]['coordinates'][0] =  Wake_Geometry[IWR]['coordinates'][0] + UWvx[IWR,0] * Dt
                Wake_Geometry[IWR]['coordinates'][2] =  Wake_Geometry[IWR]['coordinates'][2] + UWvx[IWR,1] * Dt
                
    # Tangential Velocity, pressure and loads
        for icp in range(NELEMENTS):
            Qt_upper[It,icp] = Uinf + 0.5 * Gamma[icp,It]
            Qt_lower[It,icp] = Uinf - 0.5 * Gamma[icp,It]
            DGamma[It,icp] = (Gamma[icp,It]-Gamma[icp,It-1])/Dt
            delta_P[It,icp] = rho*(Uinf * Gamma[icp,It] + np.sum(DGamma[It,:icp],axis=0))
            Drag_panel[It,icp] = WWake[It,icp]*Gamma[icp,It] + np.sum(DGamma[It,:icp],axis=0) * SN
            
        Lift[It] = np.sum(delta_P[It,:])*CS
        Drag[It] = rho * np.sum(Drag_panel[It,:],axis=0) 
        
        CL[It] = Lift[It]/(0.5 * rho * Uinf**2)
        CD[It] = Drag[It]/(0.5 * rho * Uinf**2)

        # Plot velocity field from vorticity
        Vx_list = []
        Vz_list = []
        P_list = []
        x_list = []
        z_list = []

        gammas = np.array([vortex['strenght'] for vortex in Airfoil_geometry[0]['vortexpoints']])
        vortex_x_coordinates = np.array([vortex['coordinates'][0] for vortex in Airfoil_geometry[It]['vortexpoints']] + [Wake_Vort_X])
        vortex_z_coordinates = np.array([vortex['coordinates'][1] for vortex in Airfoil_geometry[It]['vortexpoints']] + [Wake_Vort_Z])

        n_chords = 5
        c = vortex_x_coordinates[-1] - vortex_x_coordinates[0]
        x_coordinates = np.linspace(-n_chords * c + vortex_x_coordinates[0], n_chords * c + vortex_x_coordinates[0], 20)
        z_coordinates = np.linspace(-5, 5, 30)

        for x, z in itertools.product(x_coordinates, z_coordinates):
            Vx_vort, Vz_vort = vorticity_to_velocity_field(
                x,
                z,
                gammas=Gamma[:, 0],
                vortex_x_coordinates=vortex_x_coordinates,
                vortex_z_coordinates=vortex_z_coordinates
            )

            Vx = Vx_vort + Uinf
            Vz = Vz_vort

            P = 101325 - 1.225/2*(Vx**2 + Vz**2)

            Vx_list.append(Vx)
            Vz_list.append(Vz)
            P_list.append(P)
            x_list.append(x)
            z_list.append(z)

        df = pd.DataFrame({'x': x_list, 'z': z_list, 'Vx': Vx_list, 'Vz': Vz_list, 'P': P_list})

        df.sort_values(by=['x', 'z'], inplace=True)

        X, Z = np.meshgrid(x_coordinates, z_coordinates)
        # VX = df.pivot('z', columns='x', values='Vx')
        # VZ = df.pivot('z', columns='x', values='Vz')

        P = df.pivot('z', columns='x', values='P')

        if It % 20 == 0:  # plot every 20th timestep
            plt.quiver(x_list, z_list, Vx_list, Vz_list, color='g')  # plot velocity field
            # plt.streamplot(X, Z, VX, VZ, color='g')  # plot velocity field
            plt.plot([0, 5], [0, -5*SN], linewidth=3)  # plot thick line for airfoil
            plt.show()

            fig, ax = plt.subplots(1,1)
            contour = ax.contourf(X, Z, P)
            cbar = fig.colorbar(contour)
            cbar.ax.set_ylabel('Total pressure [Pa]')

            ax.plot([0, 5], [0, -5*SN], linewidth=3, color='r')  # plot thick line for airfoil
            plt.legend()
            # plt.imshow(P, origin='lower', interpolation='bilinear')
            plt.show()
            # plt.close()
    # define properties of the system
    dt=Dt
    time =np.arange(0,1.0*Nsteps,dt)

    
    # properties of the airfoil
    chord=Flat_blade_chord # chord of the airfoil
    dCn_dalpha=2.0*np.pi # lift slope
    alpha0=0.0*np.pi # alpha for which normal load is zero in steady flow
    
    # pitching motion of the airfoil
    k=.1 # reduced frequency of the pitching motion
    omega=k*2/chord*Uinf # frequency of the piching motion
    Amplitude_alpha=alpha # amplitude of the pitching motion
    alpha_t0=0/180*np.pi # alpha at time=0
    alpha=Amplitude_alpha*np.sin(omega*time)+alpha_t0 # calculate alpha
    dalpha_dt=np.gradient(alpha,time) # calculate the time derivative of alpha
    
    
    # plunge motion of the airfoil
    k_plg=.0 # reduced frequency of the plunge motion
    omega_plg=k_plg*2/chord*Uinf # frequency of the plunge motion
    Amplitude_plg=.3 # amplitude of the plunge motion
    hplg=Amplitude_plg*np.sin(omega_plg*time) #position
    dhplg_dt=np.gradient(hplg,time) # plunge velocity
    
    # define the array semi-chord time scale
    sarray = time2semichord(time, Uinf, chord)
    
    # calculate quasi-steady alpha
    # alpha0 = 0 # we define the \alpha_{0}, zero for the case of an uncambered plate/airfoil
    alphaqs = alpha + dalpha_dt*(chord/2)/Uinf - dhplg_dt/Uinf
    dalphaqs_dt=np.gradient(alphaqs,time) # calculate the time derivative of the quasi-steady alpha
    # calculate the coefficient of normal force assuming quasi-steady flow asuming potential flow
    Cnormal_quasisteady = 2*np.pi*(alphaqs-alpha0)
    
    # we plot the effective quasi-steady angle of attack \alpha_{qs}
#%%
    # plot figure
    plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
    plt.rcParams["font.family"] = "serif" # define font
    plt.rcParams["mathtext.fontset"] = "dejavuserif"  # define font
    cmap = plt.get_cmap('BuGn')  # define colormap
    fig,ax = plt.subplots(figsize=[10,10]) # define pointers for the figure and axes
    ax.plot(alpha*180/np.pi, alphaqs*180/np.pi,color='black', linewidth=1) # plot equivalent quasi-steady angle of attack
    ax.set_xlabel(r'$\alpha (^\circ)$') # set x-label
    ax.set_ylabel(r'$\alpha_{qs} (^\circ)$') # set y-label
    # add arrows to indicate the direction of the cycle
    # parr1=ax.annotate('', xy=(17.5, 20), xytext=(10,12.5), 
    #             arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
    # parr1=ax.annotate('', xy=(10, 7.5), xytext=(17.7,15), 
    #             arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
    plt.grid() # add a grid 
    # ax.set_xlim(0,30) # define limits of the axis
    # ax.set_ylim(0,30) # define limits of the axis
    plt.tight_layout() # all elements of figure inside plot area
    plt.show() # show figure
    
    filename = './figures/alpha_quasi_steady' # define name of the figure to be saved
    fig.savefig(filename+'.svg', pad_inches = 0) # save figure
    fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
    fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure
    
    
    # we plot the quasi-steady normal-force coefficient
#%%
    # plot figure
    plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
    plt.rcParams["font.family"] = "serif" # define font
    plt.rcParams["mathtext.fontset"] = "dejavuserif"  # define font
    cmap = plt.get_cmap('BuGn')  # define colormap
    fig,ax = plt.subplots(figsize=[10,10]) # define pointers for the figure and axes
    ax.plot(alpha*180/np.pi, Cnormal_quasisteady,color='black', linewidth=1) # plot equivalent quasi-steady angle of attack
    ax.set_xlabel(r'$\alpha (^\circ)$') # set x-label
    ax.set_ylabel(r'$Cn_{qs} $') # set y-label
    # add arrows to indicate the direction of the cycle
    # parr1=ax.annotate('', xy=(-2..05, -0.4), xytext=(2.1,0.16), 
    #             arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
    # parr1=ax.annotate('', xy=(1, 0.2), xytext=(-3,-0.25), 
    #             arrowprops=dict(color='black', shrink=0.05, width=.5, headwidth=3,headlength=4, linewidth=.2))
    plt.grid() # add a grid 
    # ax.set_xlim(0,30) # define limits of the axis
    # ax.set_ylim(0,3) # define limits of the axis
    plt.tight_layout() # all elements of figure inside plot area
    plt.show() # show figure
    
    filename = './figures/Cnormal_quasisteady' # define name of the figure to be saved
    # fig.savefig(filename+'.svg', pad_inches = 0) # save figure
    fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
    fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure
    
    #%%
    # define arrays for X,Y and alpha_equivalent
    Xarray=np.zeros(np.shape(time))
    Yarray=np.zeros(np.shape(time))
    
    # define the array of alpha_equivalent
    alpha_equivalent=np.zeros(np.shape(time))
    alpha_equivalent[0]=alphaqs[0]
    
    # march solution in time for alpha_E
    for i,val in enumerate(time[:-1]):
        Xarray[i+1],Yarray[i+1]=duhamel_approx(Xarray[i],Yarray[i],sarray[i+1]-sarray[i],alphaqs[i+1]-alphaqs[i])
    
    alpha_equivalent=alphaqs-Xarray-Yarray
    
    
    # plot solutions of test of duhamel_approx

    # plot figure
    plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
    plt.rcParams["font.family"] = "serif" # define font
    plt.rcParams["mathtext.fontset"] = "dejavuserif" # define font
    cmap = plt.get_cmap('BuGn') # define colormap
    fig,ax = plt.subplots(figsize=[10,10]) # define pointers for figure and axes
    
    
    #we will only plot the last cycle
    Ncycles = np.floor(time[-1]*omega/(2*np.pi)) # determine number of cycles
    n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
    i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
    i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
    i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees
    
    # plot last cycle of the simulation, steady, quasi-steady and unsteady equivalent angle of attack 
    ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1], Uinf, chord), alpha[i1:i3]*180/np.pi,color='blue',linestyle='--', label=r'$\alpha$')
    ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1], Uinf, chord), alphaqs[i1:i3]*180/np.pi,color='red',linestyle='-.', label=r'$\alpha_{qs}$')
    ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1], Uinf, chord), alpha_equivalent[i1:i3]*180/np.pi,color='green',linestyle='-', label=r'$\alpha_{eq}$')
    ax.set_xlabel('s semichords') # set x-label
    ax.set_ylabel(r'$(^\circ)$') # set y-label
    ax.set_xlim(0,2) # define limits of the axis
    # ax.set_ylim(0,30) # define limits of the axis
    ax.grid() # add grid
    ax.legend(loc='lower left')
    plt.tight_layout() # all elements of figure inside plot area
    plt.show() # show figure
    
    
    
    filename = './figures/comparison_alpha_st_qs_circ' # define name of the figure to be saved
    # fig.savefig(filename+'.svg', pad_inches = 0) # save figure
    fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
    fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure
    
    
    #%%
    # plot figure
    plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    cmap = plt.get_cmap('BuGn')
    fig,ax = plt.subplots(figsize=[10,10])
    
    #we will only plot the last cycle
    Ncycles = np.floor(time[-1]*omega/(2*np.pi))
    n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
    i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
    i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
    i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees
    
    
    ax.plot(alpha[i1:i3]*180/np.pi, alphaqs[i1:i3]*180/np.pi,color='blue',linestyle='--', label=r'$\alpha_{qs}$')
    ax.plot(alpha[i1:i3]*180/np.pi, alpha_equivalent[i1:i3]*180/np.pi,color='red',linestyle='dashdot', label=r'$\alpha_{eq}$')
    
    
    # we will plot arrows to see the direction of the cycle
    scale_arrow=3 # scale od arrow
    dx = (alpha[i1]-alpha[i1-1]) # dx of arrow
    dy = (alphaqs[i1]-alphaqs[i1-1])  # dy of arrow   
    ax.arrow(alpha[i1]*180/np.pi, alphaqs[i1]*180/np.pi, 
                 scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
                  color='blue', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle
    dx = (alpha[i2]-alpha[i2-1]) # dx of arrow
    dy = (alphaqs[i2]-alphaqs[i2-1])  # dy of arrow   
    ax.arrow(alpha[i2]*180/np.pi, alphaqs[i2]*180/np.pi, 
                 scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
                  color='blue', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle
    
    
    dx = (alpha[i1]-alpha[i1-1]) # dx of arrow
    dy = (alpha_equivalent[i1]-alpha_equivalent[i1-1])  # dy of arrow   
    ax.arrow(alpha[i1]*180/np.pi, alpha_equivalent[i1]*180/np.pi, 
                 scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
                  color='red', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle
    dx = (alpha[i2]-alpha[i2-1]) # dx of arrow
    dy = (alpha_equivalent[i2]-alpha_equivalent[i2-1])  # dy of arrow   
    ax.arrow(alpha[i2]*180/np.pi, alpha_equivalent[i2]*180/np.pi, 
                 scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
                  color='red', width=scale_arrow*.1, shape='left') # plot arrow at 0 degrees of cycle
    
    # ax.set_aspect(aspect=40.0)
    ax.set_xlabel(r'$\alpha (^\circ)$')
    ax.set_ylabel(r'$ (^\circ)$')
    # ax.set_xlim(0,time.max())
    ax.legend(loc='lower right')
    plt.grid()
    plt.tight_layout() # all elements of figure inside plot area
    
    plt.show()
    
    filename = './figures/comparison_cycle_alpha_qs_circ' # define name of the figure to be saved
    # fig.savefig(filename+'.svg', pad_inches = 0) # save figure
    fig.savefig(filename+'.pdf',pad_inches = 0) # save figure
    fig.savefig(filename+'.png', pad_inches = 0, dpi=300) # save figure
                
#%%

    # plot solutions of test of duhamel_approx
    
    # plot figure
    plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
    plt.rcParams["font.family"] = "serif" # define font
    plt.rcParams["mathtext.fontset"] = "dejavuserif" # define font
    cmap = plt.get_cmap('BuGn') # define colormap
    fig,ax = plt.subplots(figsize=[10,10]) # define pointers for figure and axes
    
    
    #we will only plot the last cycle
    Ncycles = np.floor(time[-1]*omega/(2*np.pi)) # determine number of cycles
    n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
    i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
    i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
    i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees
    
    # plot last cycle of the simulation, steady, quasi-steady and unsteady normal force coefficient  
    ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1], Uinf, chord), 
            circulatory_normal_force(2*np.pi,alpha[i1:i3],0),color='blue',linestyle='--', label=r'$Cn_{st}$')
    ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1], Uinf, chord), 
            circulatory_normal_force(2*np.pi,alphaqs[i1:i3],0),color='red',linestyle='-.', label=r'$Cn_{qs}$')
    ax.plot(time2semichord(n_of_cycle[i1:i3]-n_of_cycle[i1], Uinf, chord), 
            circulatory_normal_force(2*np.pi,alpha_equivalent[i1:i3],0),color='green',linestyle='-', label=r'$Cn_{c}$')
    ax.set_xlabel('s semichords') # set x-label
    ax.set_ylabel(r'$Cn$') # set y-label
    ax.set_xlim(0,2) # define limits of the axis
    # ax.set_ylim(0,3) # define limits of the axis
    ax.grid() # add grid
    ax.legend(loc='lower left')
    plt.tight_layout() # all elements of figure inside plot area
    plt.show() # show figure
    
    
    
    filename = './figures/comparison_Cn_st_qs_circ' # define name of the figure to be saved
    # fig.savefig(filename+'.svg', pad_inches = .0) # save figure
    fig.savefig(filename+'.pdf',pad_inches = .0) # save figure
    fig.savefig(filename+'.png', pad_inches = .0, dpi=300) # save figure
    
    #%%
    # plot figure
    plt.rcParams.update({'font.size': 14}) #, 'figure.dpi':150, 'savefig.dpi':150})
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.fontset"] = "dejavuserif"
    cmap = plt.get_cmap('BuGn')
    fig,ax = plt.subplots(figsize=[10,10])
    
    #we will only plot the last cycle
    Ncycles = np.floor(time[-1]*omega/(2*np.pi))
    n_of_cycle = time*omega/(2*np.pi) # calculate the phase of the different points of the cycle
    i1=np.argmin(np.abs(n_of_cycle-(Ncycles-1))) # index of start of cycle plotted
    i2=np.argmin(np.abs(n_of_cycle-(Ncycles-.5))) # index of 180 degrees
    i3=np.argmin(np.abs(n_of_cycle-(Ncycles))) # index of 360 degrees
    
    
    ax.plot(alpha[i1:i3]*180/np.pi, 
            circulatory_normal_force(2*np.pi,alphaqs[i1:i3],0),color='blue',linestyle='--', label=r'$Cn_{qs}$')
    ax.plot(alpha[i1:i3]*180/np.pi, 
            circulatory_normal_force(2*np.pi,alpha_equivalent[i1:i3],0),color='red',linestyle='dashdot', label=r'$Cn_{c}$')
    
    
    # we will plot arrows to see the direction of the cycle
    scale_arrow=3 # scale od arrow
    dx = (alpha[i1]-alpha[i1-1])*180/np.pi # dx of arrow
    dy = circulatory_normal_force(2*np.pi,(alphaqs[i1]-alphaqs[i1-1]),0)  # dy of arrow   
    ax.arrow(alpha[i1]*180/np.pi, circulatory_normal_force(2*np.pi,alphaqs[i1],0), 
                 scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
                  color='blue', width=scale_arrow*.02,shape='left') # plot arrow at 0 degrees of cycle
    dx = (alpha[i2]-alpha[i2-1])*180/np.pi # dx of arrow
    dy = circulatory_normal_force(2*np.pi,(alphaqs[i2]-alphaqs[i2-1]),0)  # dy of arrow   
    ax.arrow(alpha[i2]*180/np.pi, circulatory_normal_force(2*np.pi,alphaqs[i2],0), 
                 scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
                  color='blue', width=scale_arrow*.02,shape='left') # plot arrow at 0 degrees of cycle
    
    
    dx = (alpha[i1]-alpha[i1-1])*180/np.pi # dx of arrow
    dy = circulatory_normal_force(2*np.pi,(alpha_equivalent[i1]-alpha_equivalent[i1-1]),0)  # dy of arrow   
    ax.arrow(alpha[i1]*180/np.pi, circulatory_normal_force(2*np.pi,alpha_equivalent[i1],0), 
                 scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
                  color='red', width=scale_arrow*.02, shape='left') # plot arrow at 0 degrees of cycle
    dx = (alpha[i2]-alpha[i2-1])*180/np.pi # dx of arrow
    dy = circulatory_normal_force(2*np.pi,(alpha_equivalent[i2]-alpha_equivalent[i2-1]),0)  # dy of arrow   
    ax.arrow(alpha[i2]*180/np.pi, circulatory_normal_force(2*np.pi,alpha_equivalent[i2],0), 
                 scale_arrow*dx/np.sqrt(dx**2+dy**2) , scale_arrow*dy/np.sqrt(dx**2+dy**2),
                  color='red', width=scale_arrow*.02, shape='left') # plot arrow at 0 degrees of cycle
    
    # ax.set_aspect(aspect=40.0)
    ax.set_xlabel(r'$\alpha (^\circ)$')
    ax.set_ylabel(r'$Cn$')
    # ax.set_xlim(0,time.max())
    # ax.set_ylim(0,3)
    ax.legend(loc='lower right')
    plt.grid()
    plt.tight_layout() # all elements of figure inside plot area
    
    plt.show()
    
    filename = './figures/comparison_cycle_Cn_qs_circ' # define name of the figure to be saved
    # fig.savefig(filename+'.svg', pad_inches = 0.) # save figure
    fig.savefig(filename+'.pdf',pad_inches = 0.) # save figure
    fig.savefig(filename+'.png', pad_inches = 0., dpi=300) # save figure
                
                
                
                
                
                
               
        
                
                    

        
        
        
        
        
        
        
    
    
    
    
    
    
    