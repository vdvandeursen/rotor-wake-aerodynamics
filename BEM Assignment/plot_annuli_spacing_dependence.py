import matplotlib.pyplot as plt

from bem import BladeElementModel
import numpy as np
import pandas as pd

if __name__ == '__main__':

    col_names = ['angle_of_attack', 'lift_coefficient', 'drag_coefficient', 'moment_coefficient']
    airfoil = pd.read_csv('DU95W180.csv', names=col_names)

    inflow_speed = 10
    tip_speed_ratio = 6
    blade_span = 50
    blade_start = 10

    annuli_range = range(5, 50, 1)

    def twist_function(r, blade_span): return 14 * (1 - r / blade_span)
    def chord_function(r, blade_span): return 3 * (1 - r / blade_span) + 1

    results = []
    for number_of_annuli in annuli_range:
        linear_sections = np.linspace(blade_start, blade_span, number_of_annuli, endpoint=False)
        cosine_sections = blade_start + (blade_span - blade_start) * (1 - np.cos(
            np.linspace(0, np.pi/2, number_of_annuli, endpoint=False)
        ))
        sine_sections = blade_start + (blade_span - blade_start) * np.sin(
            np.linspace(0, np.pi/2, number_of_annuli, endpoint=False)
        )

        spacings = {
            'Linear': linear_sections,
            'Sine': sine_sections,
            'Cosine': cosine_sections,
        }

        for key, section_starts in spacings.items():
            blade_element_model = BladeElementModel(
                blade_span=blade_span,
                blade_start=0.2 * blade_span,
                blade_pitch=2,
                blade_number=3,
                twist=twist_function(section_starts, blade_span),
                chord=chord_function(section_starts, blade_span),
                airfoil_data=airfoil,
                section_starts=section_starts
            )

            try:
                res = blade_element_model.run_performance_calculations(
                    free_stream_velocity=inflow_speed,
                    tip_speed_ratio=tip_speed_ratio,
                    rotor_yaw=0,
                    show_plots=False,
                    air_density=1.225,
                    prandtl=True
                )
            except ValueError:
                res = {'thrust_coefficient': np.nan, 'power_coefficient': np.nan}
                print(f'Not converged for {key} for {len(section_starts)} annuli')

            results.append({
                'spacing_method': key,
                'number_of_annuli': len(section_starts),
                'thrust_coefficient': res['thrust_coefficient'],
                'power_coefficient': res['power_coefficient']
            })

    results_df = pd.DataFrame(results)
    results_df.set_index(['spacing_method', 'number_of_annuli'], inplace=True)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    axs[0].set_xlabel('Number of annuli', fontsize=12)
    axs[0].set_ylabel(r'Thrust coefficient $C_T$', fontsize=12)
    axs[0].grid()
    axs[0].set_ylim(0.1, 0.50)
    axs[1].set_xlabel('Number of annuli', fontsize=12)
    axs[1].set_ylabel(r'Power coefficient $C_T$', fontsize=12)
    axs[1].grid()
    axs[1].set_ylim(0.1, 0.50)

    for spacing in spacings:
        axs[0].plot(
            results_df.loc[spacing, 'thrust_coefficient'],
            label=spacing
        )
        axs[1].plot(
            results_df.loc[spacing, 'power_coefficient'],
            label=spacing
        )

    plt.legend()
    plt.show()
