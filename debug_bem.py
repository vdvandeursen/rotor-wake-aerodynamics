from bem import BladeElementModel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def debug_induction_factor_thrust_coefficient(blade_element_model: BladeElementModel):
    a = np.arange(-.5, 1, .01)
    CTmom = blade_element_model._calculate_thrust_coefficient(a)  # CT without correction
    CTglauert = blade_element_model._calculate_thrust_coefficient(a, True)  # CT with Glauert's correction
    a2 = blade_element_model._calculate_induction_factor(CTglauert)

    fig1 = plt.figure(figsize=(12, 6))
    plt.plot(a, CTmom, 'k-', label='$C_T$')
    plt.plot(a, CTglauert, 'b--', label='$C_T$ Glauert')
    plt.plot(a, CTglauert * (1 - a), 'g--', label='$C_P$ Glauert')
    plt.xlabel('a')
    plt.ylabel(r'$C_T$ and $C_P$')
    plt.grid()
    plt.legend()

    plt.show()

def debug_prandtl_correction_tip_root(blade_element_model: BladeElementModel):
    blade_element_model.tip_speed_ratio = 6
    blade_element_model.blade_span = 1
    blade_element_model.blade_start = 0.2
    blade_element_model.blade_number = 3

    r_R = np.arange(0.2, 1, .01)
    a = np.zeros(np.shape(r_R)) + 0.3

    Prandtl, Prandtltip, Prandtlroot = blade_element_model._prandtl_tip_root_correction(r_R, a, return_all=True)

    plt.figure(figsize=(12, 6))
    plt.plot(r_R, Prandtl, 'r-', label='Prandtl')
    plt.plot(r_R, Prandtltip, 'g.', label='Prandtl tip')
    plt.plot(r_R, Prandtlroot, 'b.', label='Prandtl root')
    plt.xlabel('r/R')
    plt.legend()
    plt.show()

def debug_cl_cd(blade_element_model):
    polar_alpha = blade_element_model.airfoil_data.angle_of_attack.to_numpy()
    polar_cl = blade_element_model.airfoil_data.lift_coefficient.to_numpy()
    polar_cd = blade_element_model.airfoil_data.drag_coefficient.to_numpy()

    # plot polars of the airfoil C-alfa and Cl-Cd

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    axs[0].plot(polar_alpha, polar_cl)
    axs[0].set_xlim([-30, 30])
    axs[0].set_xlabel(r'$\alpha$')
    axs[0].set_ylabel(r'$C_l$')
    axs[0].grid()
    axs[1].plot(polar_cd, polar_cl)
    axs[1].set_xlim([0, .1])
    axs[1].set_xlabel(r'$C_d$')
    axs[1].grid()

    plt.show()

if __name__ == '__main__':
    col_names = ['angle_of_attack', 'lift_coefficient', 'drag_coefficient', 'moment_coefficient']
    airfoil = pd.read_csv('DU95W180.csv', names=col_names)

    inflow_speed = 10
    tip_speed_ratios = [6, 8, 10]
    yaw_angles = [0, 15, 30]

    def twist_function(r, blade_span): return 14*(1-r/blade_span)
    def chord_function(r, blade_span): return 3*(1-r/blade_span)+1

    blade_element_model = BladeElementModel(
        blade_span=50,
        blade_start=0.2*50,
        blade_pitch=-2,
        blade_number=3,
        twist=twist_function,
        chord=chord_function,
        airfoil_data=airfoil
    )

    # debug_induction_factor_thrust_coefficient(blade_element_model)
    # debug_prandtl_correction_tip_root(blade_element_model)
    # debug_cl_cd(blade_element_model)

    print('break')