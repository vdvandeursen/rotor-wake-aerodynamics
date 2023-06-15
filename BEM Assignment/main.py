from bem import BladeElementModel
import itertools
import pandas as pd
import numpy as np


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
    tip_speed_ratios = np.arange(6, 12, 2)
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
        # power_coefficient, thrust_coefficient,section_loads
        
        Calculations = blade_element_model.run_performance_calculations(
            free_stream_velocity=inflow_speed,
            tip_speed_ratio=tip_speed_ratio,
            rotor_yaw=yaw_angle,
            show_plots=False,
            # save_plots_dir='./plots'
            air_density=1.225,
            prandtl=True
        )
        
        power_coefficient = Calculations['power_coefficient']
        thrust_coefficient = Calculations['thrust_coefficient']
        section_loads = Calculations['blade_section_loadings']
                
        Rotor_Torque = blade_element_model.blade_span*section_loads['rotor_tangential_force'].sum()  
        Rotor_Thrust = section_loads['rotor_axial_force'].sum()  

        results.append({
            'inflow_speed': inflow_speed,
            'tip_speed_ratio': tip_speed_ratio,
            'yaw_angle': yaw_angle,
            'power_coefficient': power_coefficient,
            'thrust_coefficient': thrust_coefficient,
            'Rotor_Torque': Rotor_Torque,
            'Rotor_Thrust': Rotor_Thrust,
        })

        blade_element_model._plot_alphaVSr_R(
            r_R=blade_element_model.blade_sections['section_start']/blade_element_model.blade_span,
            Alpha=section_loads['Alpha']
        )

    performance_df = pd.DataFrame(results)

    print(performance_df)


#
# r_R = blade_element_model.blade_sections['section_start']/blade_element_model.blade_span
# blade_element_model._plot_alphaVSr_R(r_R,section_loads['Alpha'])
# blade_element_model._plot_InflowVSr_R(r_R,section_loads['Inflow_angle'])
# blade_element_model._plot_Axial_Azmth_InductionVSr_R(r_R,section_loads['axial_induction'],section_loads['tangential_induction'])
# blade_element_model._plot_Axial_Azmth_InductionVSr_R(r_R,section_loads['axial_induction'],section_loads['tangential_induction'])
# blade_element_model._plot_Ct_CqVSr_R(r_R,section_loads['section_thrust_coefficient'],section_loads['section_torque_coefficient'])
# blade_element_model._plot_RotorThrustVSTipSpeedRatio(tip_speed_ratios,performance_df['Rotor_Thrust'])
# blade_element_model._plot_RotorTorqueVSTipSpeedRatio(tip_speed_ratios,performance_df['Rotor_Torque'])
#
