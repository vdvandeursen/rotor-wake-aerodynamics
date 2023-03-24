from bem_v2 import BladeElementModel
import itertools
import pandas as pd


"""
NOTES: 

- I've changed formula for angle of attack at airfoil: included pitch
- Added air density to lift & drag equations
- Implemented root start location for prandtl root correction



"""

if __name__ == '__main__':
    col_names = ['angle_of_attack', 'lift_coefficient', 'drag_coefficient', 'moment_coefficient']
    airfoil = pd.read_csv('DU95W180.csv', names=col_names)

    inflow_speed = 10
    tip_speed_ratios =[ 6, 7, 8, 9, 10]
    yaw_angles = [0] #, 15, 30]



    def twist_function(r, blade_span): return 14*(1-r/blade_span)
    def chord_function(r, blade_span): return 3*(1-r/blade_span)+1



    blade_element_model = BladeElementModel(
        blade_span=50,
        blade_start=0.2*50,
        blade_pitch=2,
        blade_number=3,
        twist=twist_function,
        chord=chord_function,
        airfoil_data=airfoil
    )

    results = []

    # for every combination of TSR and yaw angle: run performance calculations and append to results list
    for tip_speed_ratio, yaw_angle in itertools.product(tip_speed_ratios, yaw_angles):
        power_coefficient, thrust_coefficient,section_loads = blade_element_model.run_performance_calculations(
            free_stream_velocity=inflow_speed,
            tip_speed_ratio=tip_speed_ratio,
            rotor_yaw=yaw_angle,
            show_plots=False,
            # save_plots_dir='./plots'
            air_density=1.225
        )
        
        Rotor_Torque = blade_element_model.blade_span*section_loads['rotor_tangential_force'].sum()  
        Rotor_Thrust = section_loads['rotor_axial_force'].sum()  

        results.append({
            'inflow_speed': inflow_speed,
            'tip_speed_ratio': tip_speed_ratio,
            'yaw_angle': yaw_angle,
            'power_coefficient': power_coefficient,
            'thrust_coefficient': thrust_coefficient,
            'Rotor_Torque':Rotor_Torque,
            'Rotor_Thrust':Rotor_Thrust,
        })

    performance_df = pd.DataFrame(results)

    print(performance_df)

    results = []

    slopes=[2.9+0.01*i for i in range(5)]
    tip_chords=[0.95+0.005*i for i in range(5)]

    for slope,tip_chord in itertools.product(slopes,tip_chords):
       def chord_function(r, blade_span): return slope*(1-r/blade_span)+tip_chord

       blade_element_model = BladeElementModel(
       blade_span=50,
       blade_start=0.2*50,
       blade_pitch=2,
       blade_number=3,
       twist=twist_function,
       chord=chord_function,
       airfoil_data=airfoil
       )
    
    
       power_coefficient, thrust_coefficient , section_loads = blade_element_model.run_performance_calculations(
       free_stream_velocity=inflow_speed,
       tip_speed_ratio=tip_speed_ratio,
       rotor_yaw=yaw_angle,
       show_plots=True,
       # save_plots_dir='./plots'
       air_density=1.225
       )
    
       results.append({
       'inflow_speed': inflow_speed,
       'tip_speed_ratio': tip_speed_ratio,
       'yaw_angle': yaw_angle,
       'power_coefficient': power_coefficient,
       'thrust_coefficient': thrust_coefficient,
       'slope': slope,
       'tip_chord': tip_chord,
#       'Inflow_angle': section_loads['Inflow_angle'],
#       'Alpha': section_loads['Alpha']
       })
    
       performance_df2 = pd.DataFrame(results)
       print(performance_df2[abs(performance_df2['thrust_coefficient']-0.75)<0.05])
print('break')


#%%
import numpy as np
 
# delta_r_R = 1/499
# r_R = np.arange(0.2, 1+delta_r_R/2, delta_r_R)
r_R = blade_element_model.blade_sections['section_start']/blade_element_model.blade_span
blade_element_model._plot_alphaVSr_R(r_R,section_loads['Alpha'])
blade_element_model._plot_InflowVSr_R(r_R,section_loads['Inflow_angle'])
blade_element_model._plot_Axial_Azmth_InductionVSr_R(r_R,section_loads['axial_induction'],section_loads['tangential_induction'])
blade_element_model._plot_Axial_Azmth_InductionVSr_R(r_R,section_loads['axial_induction'],section_loads['tangential_induction'])
blade_element_model._plot_Ct_CqVSr_R(r_R,section_loads['section_thrust_coefficient'],section_loads['section_torque_coefficient'])
blade_element_model._plot_RotorThrustVSTipSpeedRatio(tip_speed_ratios,performance_df['Rotor_Thrust'])
blade_element_model._plot_RotorTorqueVSTipSpeedRatio(tip_speed_ratios,performance_df['Rotor_Torque'])

