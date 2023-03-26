from bem import BladeElementModel
import pandas as pd
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

number_of_annuli = 200
results = []


def obtain_cp_for_twist_and_chord_distributions(x0):
    t0, t1, c0, c1, c2 = x0

    section_starts = np.linspace(0.2, 1, number_of_annuli, endpoint=False)

    # twist_distribution = t0 * section_starts**2 + t1 * section_starts - (t0 + t1)  # constrain twist = 0 at tip
    # chord_distribution = c0 * (c1 - section_starts)**2 + c1 * section_starts + c2

    twist_distribution = t0 * (1 - section_starts**t1)  # this constrains twist = 0 at tip
    chord_distribution = c0 * (c1 - section_starts**c2)

    blade_element_model = BladeElementModel(
        blade_span=50,
        blade_start=0.2 * 50,
        blade_pitch=0,
        blade_number=3,
        twist=twist_distribution,
        chord=chord_distribution,
        airfoil_data=airfoil
    )

    try:
        result = blade_element_model.run_performance_calculations(
            free_stream_velocity=10,
            tip_speed_ratio=8,
            rotor_yaw=0,
            show_plots=False,
            required_thrust_coefficient=0.75,
            air_density=1.225
        )
        results.append(result)

        pd.DataFrame(results).to_csv('optimization_results.csv', index=False)

        n = len(results)

        if (n - 1) % 10 == 0:
            print(f'Trying: {x0}')
            print(f' - obtained Cp = {result["power_coefficient"]} with blade pitch = {result["blade_pitch"]} deg')

            ax1.plot(section_starts, twist_distribution + result['blade_pitch'],
                     label=f'n={n}, Cp={result["power_coefficient"]:.2f}')
            ax2.plot(section_starts, chord_distribution, label=f'n={n}, Cp={result["power_coefficient"]:.2f}')

            plt.legend()
            plt.tight_layout()
            plt.savefig(f'convergence_result_{n}.png')

        return -1 * result['power_coefficient']

    except RuntimeError:
        print(f'No convergence found')
        return 0


if __name__ == '__main__':
    col_names = ['angle_of_attack', 'lift_coefficient', 'drag_coefficient', 'moment_coefficient']
    airfoil = pd.read_csv('DU95W180.csv', names=col_names)

    fig, (ax1, ax2) = plt.subplots(2, figsize=(16,9))
    fig.suptitle('Twist and chord distribution optimization')

    ax1.set_ylabel('Twist [deg]')
    ax2.set_ylabel('Chord [m]')
    ax2.set_xlabel('Spanwise location r/R [-]')
    ax1.set_xlim(0, 1)
    ax2.set_xlim(0, 1)

    np.seterr('ignore')

    # Optimization:
    # Constrain twist and chord distributions as 2nd order polynomials with initial guess the default design.

    # t0, t1, = 0, -14,
    # c0, c1, c2 = 0, -3, 4

    # These are the values that the current design uses.
    t0, t1 = 14, 1
    c0, c1, c2 = 3, 4/3, 1

    # These result from optimization
    # t0, t1 = 13.98, 0.803
    # c0, c1, c2 = 4.19, 1.50, 0.833

    x0 = np.array([t0, t1, c0, c1, c2])

    # bounds = [
    #     (None, None),
    #     (None, 0),
    #     (None, None),
    #     (None, 0),
    #     (0, None)
    # ]

    # coefficients = opt.minimize(obtain_cp_for_twist_and_chord_distributions, x0, bounds=bounds)
    coefficients = opt.minimize(obtain_cp_for_twist_and_chord_distributions, x0)

    print('break')


