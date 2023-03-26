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
    if True:
        results = []

        #exponents=[0.005,0.02,0.05,0.1,0.2,0.6,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.5,1.6,1.8,2,4,8,16,32,64,128]
        #As=[1+i*0.05 for i in range(32)]
        #eps=0.001
        #ct_tres=0.001
        #refer to save_results.p (pickle'd performance_df2 with params above) for choice of parameters
        # As=[1.15+i*0.01 for i in range(15)]
        # exponents=[0.65+i*0.025 for i in range(18)]
        # eps=0.0001
        # ct_tres=0.0001
        #refer to save_results_2.p (pickle'd performance_df2 with params above) for choice of parameters
        As=[1.2]
        exponents=[0.725]
        eps=0.0001
        ct_tres=0.0001
        for A,exponent in itertools.product(As,exponents):
            multi=4.185976
            def chord_function(r, blade_span): return multi*(A-(r/blade_span)**(exponent))

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
            while abs(thrust_coefficient-0.75)>ct_tres:
                multi=multi+eps
                def chord_function(r, blade_span): return multi*(A-(r/blade_span)**(exponent))

                blade_element_model = BladeElementModel(
                blade_span=50,
                blade_start=0.2*50,
                blade_pitch=2,
                blade_number=3,
                twist=twist_function,
                chord=chord_function,
                airfoil_data=airfoil
                )

                power_coefficient_new, thrust_coefficient_new , section_loads_new = blade_element_model.run_performance_calculations(
                free_stream_velocity=inflow_speed,
                tip_speed_ratio=tip_speed_ratio,
                rotor_yaw=yaw_angle,
                show_plots=True,
                # save_plots_dir='./plots'
                air_density=1.225
                )

                dct=(thrust_coefficient_new-thrust_coefficient)/eps
                multi=multi-eps+(0.75-thrust_coefficient)/dct

                power_coefficient, thrust_coefficient, section_loads = blade_element_model.run_performance_calculations(
                free_stream_velocity=inflow_speed,
                tip_speed_ratio=tip_speed_ratio,
                rotor_yaw=yaw_angle,
                show_plots=True,
                # save_plots_dir='./plots'
                air_density=1.225
                )
            print(power_coefficient)


            
            results.append({
            'inflow_speed': inflow_speed,
            'tip_speed_ratio': tip_speed_ratio,
            'yaw_angle': yaw_angle,
            'power_coefficient': power_coefficient,
            'thrust_coefficient': thrust_coefficient,
            'A': A,
            'exp': exponent,
            'multi': multi
        #       'Inflow_angle': section_loads['Inflow_angle'],
        #       'Alpha': section_loads['Alpha']
            })
        
        performance_df2 = pd.DataFrame(results)
        print(performance_df2)
        print(max(performance_df2['power_coefficient']))
        print(performance_df2[performance_df2['power_coefficient']>0.99*max(performance_df2['power_coefficient'])])
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

