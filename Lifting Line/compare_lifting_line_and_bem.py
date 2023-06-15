import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifting_line_model_functions import create_rotor_geometry, solve_lifting_line_system_matrix_approach,changewakeconvectionspeedplot
from bem import BladeElementModel


"""This script creates the plots that compare the Lifting Line model and BEM model"""

#### INITIATE LIFTING LINE MODEL ####

# Generate normal wake geometry (control wake length and the element number and convection of the wake)
NELEMENTS = 10
Nrotations = 25
Rotation_resoltion = 20
aw = 0.0

s_Array = np.linspace(0, np.pi, NELEMENTS)

s_cos_array = []
for i in range(len(s_Array)):
    s_cos_array.append((-0.8/2.0 * (np.cos(s_Array[i])-1)+0.2) * 50.0 )

maxradius = max(s_cos_array)
theta_Array = np.linspace(0, Nrotations*2*np.pi, Rotation_resoltion)

TSR = 6
Uinf = 1
nblades = 1

rotor_wake_system = create_rotor_geometry(s_cos_array, maxradius, TSR/(1-aw), Uinf, theta_Array,nblades)

wind = [Uinf,0 , 0]
Omega = TSR* Uinf/maxradius
AirDen = 1.225
lift_solution = solve_lifting_line_system_matrix_approach(rotor_wake_system, wind, Omega,maxradius,nblades);

nonDim_Force = 0.5* AirDen* Uinf**2 *maxradius
nonDim_Gamma = Uinf**2*np.pi/(nblades*Omega)


#### INITIATE BEM MODEL ####
col_names = ['angle_of_attack', 'lift_coefficient', 'drag_coefficient', 'moment_coefficient']
airfoil = pd.read_csv('DU95W180.csv', names=col_names, header=0)

blade_element_model = BladeElementModel(
    blade_span=50,
    blade_start=0.2 * 50,
    blade_pitch=2,
    blade_number=3,
    twist=lambda r, R: 14 * (1 - r/R),
    chord=lambda r, R: 3 * (1 - r/R) + 1,
    airfoil_data=airfoil
)

bem_results = blade_element_model.run_performance_calculations(
        free_stream_velocity=Uinf,
        tip_speed_ratio=TSR,
        rotor_yaw=0,
        show_plots=False,
        air_density=AirDen,
        prandtl=True
    )

if __name__ == "__main__":
    # 1. Plot angle of attack vs spanwise location.
    fig = plt.figure(figsize=(9, 6))

    plt.plot(
        np.linspace(0.2, 1, 100),
        bem_results['blade_section_loadings']['Alpha'],
        label='BEM'
    )

    plt.plot(
        lift_solution['r_R'],
        lift_solution['alpha'],
        label='Lifting line'
    )

    plt.xlabel('Blade spanwise location [r/R]', fontsize=12)
    plt.ylabel('Angle of attack (deg)', fontsize=12)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.grid(axis='both')
    plt.savefig('BEM_vs_Lifting_Line_aoa.png', dpi=600)
    plt.close(fig)

    # 2. Plot induction factor vs spanwise loc
    fig = plt.figure(figsize=(9, 6))

    plt.plot(
        np.linspace(0.2, 1, 100),
        bem_results['blade_section_loadings']['axial_induction'],
        label='BEM'
    )

    plt.plot(
        lift_solution['r_R'],
        lift_solution['a'],
        label='Lifting line'
    )

    plt.xlabel('Blade spanwise location [r/R]', fontsize=12)
    plt.ylabel('Axial induction factor [-]', fontsize=12)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.grid(axis='both')
    plt.savefig('BEM_vs_Lifting_Line_induction.png', dpi=600)
    plt.close(fig)

    # 3. Plot Ct vs r.R
    fig = plt.figure(figsize=(9, 6))

    plt.plot(
        np.linspace(0.2, 1, 100),
        bem_results['blade_section_loadings']['section_thrust_coefficient'],
        label='BEM'
    )

    plt.plot(
        lift_solution['r_R'],
        lift_solution['CP_array'],
        label='Lifting line'
    )

    plt.xlabel('Blade spanwise location [r/R]', fontsize=12)
    plt.ylabel('Thrust coefficient [-]', fontsize=12)
    plt.legend(fontsize=10)
    plt.tight_layout()
    plt.grid(axis='both')
    plt.savefig('BEM_vs_Lifting_Line_thrust_coefficient.png', dpi=600)
    plt.close(fig)
