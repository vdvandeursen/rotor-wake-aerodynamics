from bem import BladeElementModel
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
    tip_speed_ratios =[8]# [6, 8, 10]
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
        power_coefficient, thrust_coefficient = blade_element_model.run_performance_calculations(
            free_stream_velocity=inflow_speed,
            tip_speed_ratio=tip_speed_ratio,
            rotor_yaw=yaw_angle,
            show_plots=False,
            # save_plots_dir='./plots'
            air_density=1.225
        )

        results.append({
            'inflow_speed': inflow_speed,
            'tip_speed_ratio': tip_speed_ratio,
            'yaw_angle': yaw_angle,
            'power_coefficient': power_coefficient,
            'thrust_coefficient': thrust_coefficient
        })

    performance_df = pd.DataFrame(results)

    print(performance_df)

    results = []

    slopes=[2.9+0.01*i for i in range(20)]
    tip_chords=[0.95+0.005*i for i in range(20)]

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


        power_coefficient, thrust_coefficient = blade_element_model.run_performance_calculations(
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
            'tip_chord': tip_chord
        })

    performance_df2 = pd.DataFrame(results)
    print(performance_df2[abs(performance_df2['thrust_coefficient']-0.75)<0.05])
print('break')
