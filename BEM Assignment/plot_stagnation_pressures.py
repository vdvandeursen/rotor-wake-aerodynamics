from bem import BladeElementModel
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

col_names = ['angle_of_attack', 'lift_coefficient', 'drag_coefficient', 'moment_coefficient']
airfoil = pd.read_csv('DU95W180.csv', names=col_names)


def twist_function(r, blade_span): return 14 * (1 - r / blade_span)
def chord_function(r, blade_span): return 3 * (1 - r / blade_span) + 1


U1 = 10
p1 = 101325
rho = 1.225

stag_p_1 = []
stag_p_2 = []
stag_p_3 = []
stag_p_4 = []
spans = []

for span in range(25, 100, 10):

    Ar = np.pi * span ** 2
    try:
        blade_element_model = BladeElementModel(
            blade_span=span,
            blade_start=0.2*span,
            blade_pitch=2,
            blade_number=3,
            twist=twist_function,
            chord=chord_function,
            airfoil_data=airfoil
        )

        res = blade_element_model.run_performance_calculations(
            free_stream_velocity=U1,
            tip_speed_ratio=6,
            rotor_yaw=0,
            show_plots=False,
            # save_plots_dir='./plots'
            air_density=1.225,
            prandtl=True,
        )

        Ft = res['thrust_coefficient'] * 0.5 * rho * U1**2 * Ar

        stag_p_1.append(p1 + 0.5 * rho * U1 ** 2)
        U4 = np.sqrt(U1**2 - 2/rho * Ft / Ar)

        p4 = p1
        stag_p_4.append(p4 + 0.5 * rho * U4 ** 2)

        Ur = 0.5 * (U1 + U4)

        p3 = p4 + rho * U4 ** 2 - rho*Ur ** 2
        stag_p_3.append(p3 + 0.5 * rho * Ur ** 2)

        p2 = Ft/Ar + p3
        stag_p_2.append(p2 + 0.5 * rho * Ur ** 2)

        spans.append(span)
    except ValueError:
        print(f'Error for {span}')
        continue

plt.figure()
plt.plot(spans, stag_p_1, label='Upstream inf')
plt.plot(spans, stag_p_2, label='Upstream rotor')
plt.plot(spans, stag_p_3, label='Downstream rotor')
plt.plot(spans, stag_p_4, label='Downstream inf')
plt.legend()
plt.xlabel('Rotor radius [m]', fontsize=12)
plt.ylabel('Stagnation pressure [Pa]', fontsize=12)

plt.show()
