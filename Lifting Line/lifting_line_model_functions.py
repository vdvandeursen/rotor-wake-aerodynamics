# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 16:44:57 2023

@author: MFekry
"""
import numpy as np
import pandas as pd



def geoBlade(r_R):
  pitch = 2
  chord = 3*(1-r_R)+1
  twist = -14*(1-r_R)
  result =[chord , twist + pitch]
  return result
  
def update_Gamma_sinle_ring(ring,GammaNew,WeightNew):
    for i in range(len(ring['filaments'])):
        ring['filaments'][i]['Gamma'] = ring['filaments'][i]['Gamma']*(1-WeightNew) + WeightNew * GammaNew
    return(ring);

def velocity_induced_single_ring(ring,controlpoint):
        tempvel1=[0 ,0, 0];
        velind=[0 ,0, 0];
        CORE = 0.00001;
        for i in range(len(ring['filaments'])):
            GAMMA = ring['filaments'][i]['Gamma'];
            XV1 = [ ring['filaments'][i]['x1'] , ring['filaments'][i]['y1']   ,  ring['filaments'][i]['z1'] ];
            XV2 = [ ring['filaments'][i]['x2'] , ring['filaments'][i]['y2']   ,  ring['filaments'][i]['z2'] ];
            # XVP = controlpoint;
            tempvel1 =velocity_3D_from_vortex_filament(GAMMA,XV1, XV2, controlpoint ,CORE);
            velind[0]+=tempvel1[0];
            velind[1]+=tempvel1[1];
            velind[2]+=tempvel1[2];
            # // console.log(velind)
        return velind;
    
def polarAirfoil(alpha):
    airfoil = 'DU95W180.csv'
    data1=pd.read_csv(airfoil, header=0, names = ["alfa", "cl", "cd", "cm"],  sep=',')
    polar_alpha = data1['alfa'][:]
    polar_cl = data1['cl'][:]
    polar_cd = data1['cd'][:]
    cl = np.interp(alpha,polar_alpha,polar_cl);
    cd = np.interp(alpha,polar_alpha,polar_cd);
    result =[cl , cd];
    return result;
 

def calculateCT_CProtor(a,aline,Fnorm,Ftan, Uinf, r_Rarray, Omega, Radius, NBlades):
    # calculates the performance of the rotor, retunring CT and CPs
    r_R_temp=0; # radial position for evaluation
    drtemp=0; # delta r
    CProtor = 0;
    CTrotor = 0;
    CP_array = np.zeros((len(r_Rarray),1))
    CT_array = np.zeros((len(r_Rarray),1))
    for i in range(len(r_Rarray)-1):
        r_R_temp = (r_Rarray[i]+r_Rarray[i+1])/2;
        drtemp = (-r_Rarray[i]+r_Rarray[i+1]);
        CT_array[i+1] = ((drtemp*Fnorm[i]*NBlades)/(0.5*Uinf**2*np.pi*Radius*drtemp))
        CP_array[i+1] = ((drtemp*Ftan[i]*r_R_temp*Omega*NBlades)/(0.5*Uinf**3*np.pi*drtemp))
        CTrotor+=(drtemp*Fnorm[i]*NBlades)/(0.5*Uinf**2*np.pi*Radius);
        CProtor+=(drtemp*Ftan[i]*r_R_temp*Omega*NBlades)/(0.5*Uinf**3*np.pi);
    results = [CTrotor,CT_array, CProtor,CP_array];
    return results;
     
def changewakeconvectionspeedplot(aw,Uinf,TSR,nblades):
    s_Array = np.linspace(0, np.pi, 10)
    for i in range(len(s_Array)):
        s_Array[i]= (-0.8/2.0 * (np.cos(s_Array[i])-1)+0.2) * 50.0;
        maxradius = max(s_Array);
        theta_Array = np.linspace(0, 50*np.pi, 20)
        rotor_wake_system = create_rotor_geometry(s_Array, maxradius, TSR/(1-aw), Uinf, theta_Array,nblades);
    return rotor_wake_system





def loadBladeElement(Vnorm, Vtan, r_R):
  Vmag2 = (pow(Vnorm,2) + pow(Vtan,2)  );
  InflowAngle = np.arctan(Vnorm/Vtan);
  temp= geoBlade(r_R);
  chord = temp[0];
  twist = temp[1];
  alpha = twist + InflowAngle*180/np.pi;
  temp = polarAirfoil(alpha);
  cl = temp[0];
  cd = 0*temp[1];

  Lift = 0.5*Vmag2*cl*chord;
  Drag = 0.5*Vmag2*cd*chord;
  Fnorm = Lift*np.cos(InflowAngle)+Drag*np.sin(InflowAngle);
  Ftan = Lift*np.sin(InflowAngle)-Drag*np.cos(InflowAngle);
  Gamma = 0.5*np.sqrt(Vmag2)*cl*chord ; 
  result = [Fnorm , Ftan, Gamma, alpha, InflowAngle];
  return result;





def create_rotor_geometry(span_array, radius, tipspeedratio, Uinf, theta_array, nblades):
    
  filaments = []
  ring = []
  controlpoints = []
  bladepanels = []

  for krot in range(nblades): 

    angle_rotation = 2*np.pi/nblades*krot
    cosrot = np.cos(angle_rotation)
    sinrot = np.sin(angle_rotation)

    for i in range(len(span_array)-1): 
      r = (span_array[i]+span_array[i+1])/2;
      geodef = geoBlade(r/radius);
      angle = geodef[1]*np.pi/180;
      # define controlpoints
      temp1 = {'coordinates': [ 0 ,  r  , 0 ] , 'chord': geodef[0], 'normal': [ np.cos(angle),  0, -1*np.sin(angle)] , 'tangential': [-1*np.sin(angle), 0, -1*np.cos(angle)]   }
      # rotate blade to position
      temp1['coordinates'] = [ 0 ,  temp1['coordinates'][1]*cosrot -  temp1['coordinates'][2]*sinrot , temp1['coordinates'][1]*sinrot +  temp1['coordinates'][2]*cosrot ]
      temp1['normal'] = [ temp1['normal'][0],  temp1['normal'][1]*cosrot -  temp1['normal'][2]*sinrot , temp1['normal'][1]*sinrot +  temp1['normal'][2]*cosrot ]
      temp1['tangential'] = [ temp1['tangential'][0] ,  temp1['tangential'][1]*cosrot -  temp1['tangential'][2]*sinrot , temp1['tangential'][1]*sinrot +  temp1['tangential'][2]*cosrot ]

      controlpoints.append( temp1 );
      # define bound vortex filament
      temp1= {'x1':0 , 'y1':span_array[i], 'z1':0, 'x2':0 , 'y2':span_array[i+1], 'z2':0, 'Gamma': 0  }   
      # rotate filament to position
      filaments.append(temp1);
      # create trailing filaments, at x1 of bound filament
      geodef = geoBlade(span_array[i]/radius)
      angle = geodef[1]*np.pi/180
      temp1= {'x1': geodef[0]*np.sin(-1*angle), 'y1':span_array[i], 'z1': -1*geodef[0]*np.cos(angle), 'x2':0 , 'y2':span_array[i], 'z2':0, 'Gamma': 0  }  
      filaments.append(temp1);
      for j in range(len(theta_array)-1):
        xt = filaments[len(filaments)-1]['x1']
        yt = filaments[len(filaments)-1]['y1']
        zt = filaments[len(filaments)-1]['z1']
        dy = (np.cos(-theta_array[j+1])-np.cos(-theta_array[j])) * span_array[i]
        dz = (np.sin(-theta_array[j+1])-np.sin(-theta_array[j])) * span_array[i]
        dx = (theta_array[j+1]-theta_array[j])/tipspeedratio*radius

        temp1= {'x1': xt+dx, 'y1': yt+dy, 'z1': zt+dz, 'x2': xt , 'y2':yt, 'z2':zt, 'Gamma': 0  }   
        # rotate filament to position
        filaments.append(temp1);
      
      # create trailing filaments, at x2 of bound filament
      geodef = geoBlade(span_array[i+1]/radius);
      angle = geodef[1]*np.pi/180;
      temp1= {'x1':0 , 'y1':span_array[i+1], 'z1': 0, 'x2':geodef[0]*np.sin(-1*angle) , 'y2':span_array[i+1], 'z2':-1*geodef[0]*np.cos(angle), 'Gamma': 0  }  
      filaments.append(temp1)
      for j in range(len(theta_array)-1):
        xt = filaments[len(filaments)-1]['x2']
        yt = filaments[len(filaments)-1]['y2']
        zt = filaments[len(filaments)-1]['z2']
        dy = (np.cos(-theta_array[j+1])-np.cos(-theta_array[j])) * span_array[i+1]
        dz = (np.sin(-theta_array[j+1])-np.sin(-theta_array[j])) * span_array[i+1]
        dx = (theta_array[j+1]-theta_array[j])/tipspeedratio*radius

        temp1= {'x1': xt, 'y1': yt, 'z1': zt, 'x2': xt+dx , 'y2':yt+dy, 'z2':zt+dz, 'Gamma': 0  }   
        # rotate filament to position

        filaments.append(temp1);
      

      for ifil in range(len(filaments)):
        temp1=filaments[ifil]
        temp2 = [ temp1['y1']*cosrot -  temp1['z1']*sinrot , temp1['y1']*sinrot +  temp1['z1']*cosrot , temp1['y2']*cosrot -  temp1['z2']*sinrot , temp1['y2']*sinrot +  temp1['z2']*cosrot ]
        temp1['y1'] = temp2[0]
        temp1['z1'] = temp2[1]
        temp1['y2'] = temp2[2]
        temp1['z2'] = temp2[3]
        filaments[ifil] = temp1
      

      ring.append({'filaments': filaments});
      filaments = [];

      # panel of the blade section
      geodef = geoBlade(span_array[i]/radius);
      angle = geodef[1]*np.pi/180;
      geodef2 = geoBlade(span_array[i+1]/radius);
      angle2 = geodef2[1]*np.pi/180;

      temp1= {
        'p1': [-0.25*geodef[0]*np.sin(-1*angle) , span_array[i], 0.25*geodef[0]*np.cos(angle)],
        'p2': [-0.25*geodef2[0]*np.sin(-1*angle2) , span_array[i+1], 0.25*geodef2[0]*np.cos(angle2)],
        'p3': [0.75*geodef2[0]*np.sin(-1*angle2) , span_array[i+1], -0.75*geodef2[0]*np.cos(angle2)],
        'p4': [0.75*geodef[0]*np.sin(-1*angle) , span_array[i], -0.75*geodef[0]*np.cos(angle)]
      }
      temp1['p1'] = [ 0 ,  temp1['p1'][1]*cosrot -  temp1['p1'][2]*sinrot , temp1['p1'][1]*sinrot +  temp1['p1'][2]*cosrot ]
      temp1['p2'] = [ 0 ,  temp1['p2'][1]*cosrot -  temp1['p2'][2]*sinrot , temp1['p2'][1]*sinrot +  temp1['p2'][2]*cosrot ]
      temp1['p3'] = [ 0 ,  temp1['p3'][1]*cosrot -  temp1['p3'][2]*sinrot , temp1['p3'][1]*sinrot +  temp1['p3'][2]*cosrot ]
      temp1['p4'] = [ 0 ,  temp1['p4'][1]*cosrot -  temp1['p4'][2]*sinrot , temp1['p4'][1]*sinrot +  temp1['p4'][2]*cosrot ]

      bladepanels.append(temp1);
    
  return({'controlpoints': controlpoints ,  'rings': ring, 'bladepanels':bladepanels})

def solve_lifting_line_system_matrix_approach(rotor_wake_system,wind, Omega, rotorradius,nblades):
      # this codes solves a lifting line model of a horizontal axis rotor
      # as inputs, it takes
      #       rotor_wake_system: data structure that contains the geometry of the horseshoe vortex rings,
      #                         and the control points at the blade
      #       wind: unperturbed wind velocity, also known as U_infinity
      #       Omega: rotational velocity of the rotor
      #       rotorradius: the radius of the rotor

      # get controlpoints data structure
      controlpoints = rotor_wake_system['controlpoints']
      # get horseshoe vortex rings data structure
      rings = rotor_wake_system['rings']
      
      # initialize variables that we will use during the calculation
      velocity_induced =[] # velocity induced by a horse vortex ring at a control point
      GammaNew=np.zeros((len(controlpoints),1)); #// new estimate of bound circulation
      Gamma=np.zeros((len(controlpoints),1)); #// current solution of bound circulation
      for i in range(len(controlpoints)): 
          GammaNew[i]=(0) #  initialize as zeros
      MatrixU = np.zeros((len(controlpoints),len(rings))); # matrix of induction, for velocity component in x-direction
      MatrixV = np.zeros((len(controlpoints),len(rings))); # matrix of induction, for velocity component in y-direction
      MatrixW = np.zeros((len(controlpoints),len(rings))); # matrix of induction, for velocity component in z-direction
      # output variables
      a_temp = np.zeros((len(controlpoints),1)); # output vector for axial induction
      aline_temp = np.zeros((len(controlpoints),1));#   output vector for azimuthal induction
      r_R_temp = np.zeros((len(controlpoints),1));  # output vector for radial position
      Fnorm_temp = np.zeros((len(controlpoints),1)); #  output vector for axial force
      Ftan_temp = np.zeros((len(controlpoints),1));  # output vector for tangential force
      Gamma_temp = np.zeros((len(controlpoints),1)); #  output vector for circulation
      alpha_temp = np.zeros((len(controlpoints),1));
      inflow_temp = np.zeros((len(controlpoints),1));
      
      
      # the variables below are to setup the maximum number of iterations and convergence criteria
      Niterations =12000;
      errorlimit = 0.01;
      error = 1.0; #var refererror;
      ConvWeight =0.3;

      #  initalize and calculate matrices for velocity induced by horseshoe vortex rings
      # two "for cicles", each line varying wind controlpoint "icp", each column varying with
      # horseshoe vortex ring "jring"
      for icp in range(len(controlpoints)): 
          for jring in range(len(rings)):
              # set ring strenth to unity, to calculate velocity induced by horseshoe vortex ring "jring"
              # at controlpoint "icp"
              rings[jring] = update_Gamma_sinle_ring(rings[jring],1,1);
              velocity_induced = velocity_induced_single_ring(rings[jring], controlpoints[icp]['coordinates']);
              # add compnent of velocity per unit strength of circulation to induction matrix
              MatrixU[icp][jring] = velocity_induced[0];
              MatrixV[icp][jring] = velocity_induced[1];
              MatrixW[icp][jring] = velocity_induced[2];
        

      # calculate solution through an iterative process
      for kiter in range(Niterations):
          print('kiter',kiter)
          # Gamma =[]
          for ig in range(len(GammaNew)):
              # print('GammaNew:', len(GammaNew))
              # print('ig:' ,ig)
              Gamma[ig]=(GammaNew[ig]);# //update current bound circulation with new estimate
          # calculate velocity, circulation and loads at the controlpoints
          for icp in range(len(controlpoints)): 
              # print('icp:' ,icp)
              # determine radial position of the controlpoint;
              radialposition = np.sqrt(np.dot(controlpoints[icp]['coordinates'], controlpoints[icp]['coordinates']));
              u=0; v=0; w=0; #// initialize velocity
              # multiply icp line of Matrix with vector of circulation Gamma to calculate velocity at controlpoint
              for jring in range(len(rings)):
                  u = u + MatrixU[icp][jring]*Gamma[jring];# // axial component of velocity
                  v = v + MatrixV[icp][jring]*Gamma[jring]; #// y-component of velocity
                  w = w + MatrixW[icp][jring]*Gamma[jring]; #// z-component of velocity
              
              # calculate total perceived velocity
              vrot = np.cross([-Omega, 0 , 0]  , controlpoints[icp]['coordinates'] ); #// rotational velocity
              vel1 = [wind[0]+ u + vrot[0], wind[1]+ v + vrot[1] , wind[2]+ w + vrot[2]]; #// total perceived velocity at section
              # calculate azimuthal and axial velocity
              azimdir = np.cross([-1/radialposition, 0 , 0]  , controlpoints[icp]['coordinates'] ); #// rotational direction
              vazim = np.dot(azimdir , vel1); # // azimuthal direction
              vaxial =  np.dot([1, 0, 0] , vel1); #// axial velocity
              # calculate loads using blade element theory              
              temploads = loadBladeElement(vaxial, vazim, radialposition/rotorradius);
              # new point of new estimate of circulation for the blade section
              GammaNew[icp] = temploads[2];
              # update output vector
              a_temp[icp]=((-(u + vrot[0])/wind[0]));
              # print('a_temp:',np.shape(a_temp))
              aline_temp[icp]=((vazim/(radialposition*Omega)-1));
              # print('aline_temp:',np.shape(aline_temp))
              r_R_temp[icp]=(radialposition/rotorradius);
              # print('r_R_temp:',np.shape(r_R_temp))
              Fnorm_temp[icp]=(temploads[0]);
              # print('Fnorm_temp:',np.shape(Fnorm_temp))
              Ftan_temp[icp]=(temploads[1]);
              # print('Ftan_temp:',np.shape(Ftan_temp))
              Gamma_temp[icp]=(temploads[2]);
              # print('Gamma_temp:',np.shape(Gamma_temp))
              alpha_temp[icp]=(temploads[3]);
              # print('alpha_temp:',np.shape(alpha_temp))
              inflow_temp[icp]=(temploads[4]);
              # print('inflow_temp:',np.shape(inflow_temp))
              # end loop control points

          # check convergence of solution
          # print('GammaNew:',GammaNew)
          # print('Gamma:',Gamma)
          refererror =max(np.absolute(GammaNew));
          refererror =max(refererror,0.001);#  define scale of bound circulation
          error =max(np.absolute(np.subtract(np.array(GammaNew), np.array(Gamma))));# // difference betweeen iterations
          error= error/refererror;#  relative error
          # print('Error:', error)
          print('error:',error)
          if (error < errorlimit):
            #  if error smaller than limit, stop iteration cycle
            print('BREAK..............')
            break
        
      # set new estimate of bound circulation
      for ig in range(len(GammaNew)):
          GammaNew[ig] = (1-ConvWeight)*Gamma[ig] + ConvWeight*GammaNew[ig];
            
          # end iteration loop
          
          # calculate CT and CP
          
          CT_CP = calculateCT_CProtor(a_temp,aline_temp,Fnorm_temp,Ftan_temp, wind[0], r_R_temp, Omega, rotorradius, nblades)
          
          
      #  output results of converged solution
      return({'a': a_temp , 'aline': aline_temp, 'r_R': r_R_temp, 'Fnorm': Fnorm_temp, 'Ftan': Ftan_temp , 'Gamma': Gamma_temp, 'CT': CT_CP[0], 'CP': CT_CP[2], 'CT_array': CT_CP[1], 'CP_array': CT_CP[3], 'alpha': alpha_temp, 'inflow': inflow_temp});


# 3D velocity induced by a vortex filament
def velocity_3D_from_vortex_filament(GAMMA,XV1, XV2, XVP1,CORE):
    
    #  read coordinates that define the vortex filament
    X1 =XV1[0]; Y1 =XV1[1]; Z1 =XV1[2]; # // start point of vortex filament
    X2 =XV2[0]; Y2 =XV2[1]; Z2 =XV2[2]; # // end point of vortex filament
    #  read coordinates of target point where the velocity is calculated
    XP =XVP1[0]; YP =XVP1[1]; ZP =XVP1[2];
    #  calculate geometric relations for integral of the velocity induced by filament
    R1= np.sqrt(pow((XP-X1), 2) + pow((YP-Y1), 2) + pow((ZP-Z1), 2) );
    R2= np.sqrt(pow((XP-X2), 2) + pow((YP-Y2), 2) + pow((ZP-Z2), 2) );
    R1XR2_X=(YP-Y1)*(ZP-Z2)-(ZP-Z1)*(YP-Y2);
    R1XR2_Y=-(XP-X1)*(ZP-Z2)+(ZP-Z1)*(XP-X2);
    R1XR2_Z=(XP-X1)*(YP-Y2)-(YP-Y1)*(XP-X2);
    R1XR_SQR=pow(R1XR2_X, 2)+ pow(R1XR2_Y, 2)+ pow(R1XR2_Z, 2);
    R0R1 = (X2-X1)*(XP-X1)+(Y2-Y1)*(YP-Y1)+(Z2-Z1)*(ZP-Z1);
    R0R2 = (X2-X1)*(XP-X2)+(Y2-Y1)*(YP-Y2)+(Z2-Z1)*(ZP-Z2);
    #  check if target point is in the vortex filament core,
    #  and modify to solid body rotation
    if (R1XR_SQR < pow(CORE,2)):
        R1XR_SQR = pow(CORE,2);
        # GAMMA = 0;
    if (R1 < CORE):
        R1 = CORE;
        # GAMMA = 0;
    if (R2 < CORE):
        R2=CORE;
        # GAMMA = 0;
        
        # determine scalar
    K = GAMMA/4/np.pi/R1XR_SQR*(R0R1/R1 - R0R2/R2 );
    # determine the three velocity components
    U=K*R1XR2_X;
    V=K*R1XR2_Y;
    W=K*R1XR2_Z;
    # output results, vector with the three velocity components
    results = [U, V, W];
    return(results)