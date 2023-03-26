import pandas as pd
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from typing import Union, Callable

number_of_annuli = 500


# class FlowCase:
#     def __init__(self, free_stream_velocity, yaw_angle, tip_speed_ratio):
#         self.free_stream_velocity = free_stream_velocity
#         self.yaw_angle = yaw_angle
#         self.tip_speed_ratio = tip_speed_ratio


class BladeElementModel:
    def __init__(self,
                 blade_span: float,
                 blade_start: float,
                 blade_number: int,
                 blade_pitch: float,
                 airfoil_data: pd.DataFrame,
                 twist: Union[npt.ArrayLike, Callable],
                 chord: Union[npt.ArrayLike, Callable],
                 ):
        """
        Creates BladeElementModel instance. Initiates model configuration with inflow conditions and blade geometry.

        :param blade_span: blade span [m].
        :param blade_start: spanwise dimensionless location at which airfoil starts [r/R].
        :param blade_number: number of blades [-].
        :param blade_pitch: pitch angle of blade.
        :param airfoil_data: Pandas.DataFrame containing lift, drag and moment coefficients of airfoil at varying aoa.
        :param twist: specifies twist distribution. Can be a function or a numpy array (of equal length as chord array)
        :param chord: specifies chord distribution. Can be a function or a numpy array (of equal length as twist array)
        """
        self.blade_span = blade_span
        self.blade_start = blade_start
        self.blade_number = blade_number
        self.blade_pitch = blade_pitch

        self.airfoil_data = airfoil_data

        # check if twist and chord are specified as functions or as numpy arrays
        if callable(twist) and callable(chord):
            section_starts = np.linspace(blade_start, blade_span, number_of_annuli, endpoint=False)
            twist_distribution = twist(section_starts, blade_span)
            chord_distribution = chord(section_starts, blade_span)

        elif len(twist) == len(chord):
            section_starts = np.linspace(blade_start, blade_span, len(twist), endpoint=False)
            twist_distribution = twist
            chord_distribution = chord

        else:
            raise ValueError('Twist and chord must either both be specified as functions or as arrays of equal length.')

        # define blade sections as a pandas DataFrame - allows for easier calculations later on.
        self.blade_sections = pd.DataFrame(
            {
                'section_start': section_starts,
                'twist': twist_distribution,
                'chord': chord_distribution
            }
        )

        self.section_thickness = section_starts[1] - section_starts[0]  # annulus thickness

        self.free_stream_velocity = None
        self.air_density = None
        self.tip_speed_ratio = None
        self.rotor_yaw = None
        self.angular_velocity = None

    def run_performance_calculations(self,
                                     free_stream_velocity: float,
                                     tip_speed_ratio: float,
                                     rotor_yaw: float,
                                     air_density: float = 1.225,
                                     show_plots=False,
                                     save_plots_dir: str = None,
                                     prandtl_flag1 = True) -> tuple:
        """
        Perform BEM calculations with given inflow conditions.

        :param free_stream_velocity: inflow speed [m/s]
        :param air_density: inflow air density [kg/m3]
        :param tip_speed_ratio: tip speed ratio [-]
        :param rotor_yaw: rotor yaw [deg]
        :param show_plots: display performance plots
        :param save_plots_dir: if specified, saves performance plots in this directory

        :returns power_coefficient, thrust_coefficient
        """

        self.free_stream_velocity = free_stream_velocity
        self.air_density = air_density
        self.tip_speed_ratio = tip_speed_ratio
        self.rotor_yaw = rotor_yaw
        self.angular_velocity = tip_speed_ratio * free_stream_velocity / self.blade_span

        section_loadings = self.blade_sections.apply(
            lambda x: self._solve_stream_tube_for_blade_element(
                section_start=x['section_start'],
                twist=x['twist'],
                chord=x['chord'], 
                prandtl_flag = prandtl_flag1
            ),
            axis=1
        )

        dynamic_pressure = 0.5 * self.air_density * self.free_stream_velocity ** 2
        disk_surface = np.pi * self.blade_span ** 2

        rotor_thrust_coefficient = section_loadings['rotor_axial_force'].sum() / dynamic_pressure / disk_surface
        rotor_power_coefficient = (
            section_loadings['rotor_tangential_force'] *
            self.blade_sections['section_start'] *
            self.angular_velocity / (
                    dynamic_pressure * disk_surface * self.free_stream_velocity
            )
        ).sum()

        if save_plots_dir is not None:
            plt.savefig(f'{save_plots_dir}.png')

        if show_plots:
            plt.show()

        return rotor_power_coefficient, rotor_thrust_coefficient , section_loadings

    def _solve_stream_tube_for_blade_element(self, section_start, twist, chord,prandtl_flag = True):
        """
        Solve balance of momentum between blade element load and loading in the stream tube.

        :param section_start: section start location
        :param twist: section twist angle [deg]
        :param chord: section chord length [m]
        """
        stream_tube_area = np.pi * ((section_start + self.section_thickness) ** 2 - section_start ** 2)

        # initialize variables
        axial_induction_factor = 0  # axial induction
        tangential_induction_factor = 0  # tangential induction factor

        error_threshold = 0.0001  # error limit for iteration pocess, in absolute value of induction
        iteration_error = np.inf
        iteration_limit = 10e3

        result = None

        while iteration_error > error_threshold:
            # Calculate velocity and loads at blade element
            axial_velocity_rotor = self.free_stream_velocity * (1 - axial_induction_factor)
            tangential_velocity_rotor = (1 + tangential_induction_factor) * self.angular_velocity * section_start

            blade_axial_force, blade_tangential_force, gamma, inflow_angle, alpha = self._get_blade_element_loading(
                axial_velocity=axial_velocity_rotor,
                tangential_velocity=tangential_velocity_rotor,
                twist=twist,
                chord=chord
            )  # N/m

            rotor_axial_force = self.blade_number * blade_axial_force * self.section_thickness  # N
            rotor_tangential_force = self.blade_number * blade_tangential_force * self.section_thickness  # N

            # Update estimation for induction factors
            section_thrust_coefficient = rotor_axial_force / (
                    0.5 * self.air_density * stream_tube_area * self.free_stream_velocity ** 2
            )
            
            section_torque_coefficient = self.blade_span * rotor_tangential_force / (
                    0.5 * self.air_density * stream_tube_area * self.free_stream_velocity ** 2
            )

            # Update induction, accounting for Glauert's correction
            axial_induction_factor_update = self._calculate_induction_factor(section_thrust_coefficient)

            # correct new axial induction with Prandtl's correction
            prandtl_correction = self._prandtl_tip_root_correction(
                location=section_start + self.section_thickness * 0.5,  # calculate correction @ section midpoint
                axial_induction=axial_induction_factor_update
            )

            if prandtl_correction < 0.0001:
                prandtl_correction = 0.0001
                
            if prandtl_flag == False:
                prandtl_correction = 1.0
                

            axial_induction_factor_update = axial_induction_factor_update / prandtl_correction

            # for improving convergence, weigh current and previous iteration of axial induction
            axial_induction_factor_update = 0.75 * axial_induction_factor + 0.25 * axial_induction_factor_update

            # update estimate for tangential induction factor (and correcting w/ Prandtl)
            tangential_induction_factor_update = blade_tangential_force * self.blade_number / (
                2 * np.pi * self.air_density * self.free_stream_velocity * (1 - axial_induction_factor_update) *
                self.angular_velocity * 2 * section_start ** 2 * prandtl_correction
            )

            # update iteration error
            iteration_error = np.abs(axial_induction_factor_update - axial_induction_factor)

            iteration_limit -= 1

            if iteration_limit < 0:
                raise ValueError(f'Iteration threshold reached - no convergence obtained.')

            axial_induction_factor = axial_induction_factor_update
            tangential_induction_factor = tangential_induction_factor_update

            result = {
                'axial_induction': axial_induction_factor,
                'tangential_induction': tangential_induction_factor,
                'rotor_axial_force': rotor_axial_force,
                'rotor_tangential_force': rotor_tangential_force,
                'blade_axial_force': blade_axial_force,
                'blade_tangential_force': blade_tangential_force,
                "Inflow_angle": inflow_angle,
                "Alpha":alpha,
                'section_thrust_coefficient': section_thrust_coefficient,
                'section_torque_coefficient': section_torque_coefficient,
            }

        return pd.Series(result)

    def _get_blade_element_loading(self, axial_velocity, tangential_velocity, chord, twist) -> tuple:
        """
        Calculates the load for an arbitrary blade element in certain flow conditions.

        :param axial_velocity: velocity in normal direction
        :param tangential_velocity: velocity in tangential direction
        :param chord: chord length of blade section
        :param twist: twist angle of blade section

        :return axial_force, tangential_force, gamma. Forces are in [N/m]

        """
        velocity_squared = axial_velocity ** 2 + tangential_velocity ** 2
        inflow_angle = np.arctan2(axial_velocity, tangential_velocity)  # deg
        alpha = inflow_angle * 180 / np.pi - twist + self.blade_pitch  # deg

        lift_coefficient = np.interp(alpha, self.airfoil_data['angle_of_attack'], self.airfoil_data['lift_coefficient'])
        drag_coefficient = np.interp(alpha, self.airfoil_data['angle_of_attack'], self.airfoil_data['drag_coefficient'])

        lift = 0.5 * velocity_squared * lift_coefficient * chord * self.air_density  # N/m
        drag = 0.5 * velocity_squared * drag_coefficient * chord * self.air_density  # N/m

        axial_force = lift * np.cos(inflow_angle) + drag * np.sin(inflow_angle)
        tangential_force = lift * np.sin(inflow_angle) - drag * np.cos(inflow_angle)
        gamma = 0.5 * np.sqrt(velocity_squared) * lift_coefficient * chord

        return axial_force, tangential_force, gamma, inflow_angle, alpha

    def _prandtl_tip_root_correction(self, location, axial_induction, return_all=False):
        """
        Calculates Prandtl's correction based on section location and axial induction.

        :param location: location at which to calculate correction
        :param axial_induction: axial induction factor at location
        :param return_all: whether to return the root and tip corrections separately

        :return prandtl_correction

        """

        d = 2 * np.pi * (1 - axial_induction) / (
            self.blade_number * np.sqrt(self.tip_speed_ratio ** 2 + (1 - axial_induction)**2)
        )

        tip_correction = 2 / np.pi * np.arccos(
            np.exp(-np.pi * (1 - location / self.blade_span) / d)
        )

        root_correction = 2 / np.pi * np.arccos(
            np.exp(-np.pi * ((location - self.blade_start) / self.blade_span) / d)
        )

        prandtl_correction = tip_correction * root_correction

        if return_all:
            return prandtl_correction, tip_correction, root_correction

        return prandtl_correction

    @staticmethod
    def _calculate_thrust_coefficient(induction_factor, glauert_correction=False):
        """
        This function calculates the thrust coefficient as a function of induction factor 'a'

        :param glauert_correction: whether Glauert correction for heavily loaded rotors should be used.
        :param induction_factor: induction factor.
        """

        ct = 4 * induction_factor * (1 - induction_factor)

        if glauert_correction:
            ct1 = 1.816
            a1 = 1 - np.sqrt(ct1) / 2
            ct[induction_factor > a1] = ct1 - 4 * (np.sqrt(ct1) - 1) * (1 - induction_factor[induction_factor > a1])

        return ct

    @staticmethod
    def _calculate_induction_factor(ct):
        """
        This function calculates the induction factor 'a' as a function of thrust coefficient CT
        including Glauert's correction
        """
        a = np.zeros(np.shape(ct))
        ct1 = 1.816
        ct2 = 2 * np.sqrt(ct1) - ct1
        a[ct >= ct2] = 1 + (ct[ct >= ct2] - ct1) / (4 * (np.sqrt(ct1) - 1))
        a[ct < ct2] = 0.5 - 0.5 * np.sqrt(1 - ct[ct < ct2])

        return a
    
    @staticmethod
    def _plot_alphaVSr_R(r_R,Alpha):
        plt.figure(figsize=(12, 6))
        plt.title('Alpha vs r/R')
        plt.plot(r_R,Alpha)
        plt.grid()
        plt.xlabel('r/R')
        plt.ylabel('Alpha (deg)')
        plt.show()
        plt.savefig("Alpha vs r_R.jpg", bbox_inches='tight')
        plt.close('all')
        
    @staticmethod
    def _plot_InflowVSr_R(r_R,Inflow):
        plt.figure(figsize=(12, 6))
        plt.title('Inflow Angle vs r/R')
        plt.plot(r_R,Inflow)
        plt.grid()
        plt.xlabel('r/R')
        plt.ylabel('Inflow (deg)')
        plt.show()
        plt.savefig("Inflow Angle vs r_R.jpg", bbox_inches='tight')
        plt.close('all')
        
    @staticmethod
    def _plot_Axial_Azmth_InductionVSr_R(r_R,a,a_dash):
        plt.figure(figsize=(12, 6))
        plt.title('Axial Azimuthal Induction vs r/R')
        plt.plot(r_R,a,label = 'a')
        plt.plot(r_R,a_dash,label = r'$a^,$')
        plt.grid()
        plt.xlabel('r/R')
        #plt.ylabel('a')
        plt.legend()
        plt.show()
        plt.savefig("Axial Azimuthal Induction vs r_R.jpg", bbox_inches='tight')
        plt.close('all')
        
    @staticmethod
    def _plot_Ct_CqVSr_R(r_R,Ct,Cq):
        Ct[0] = 0.0
        Cq[0] = 0.0
        plt.figure(figsize=(12, 6))
        plt.title('Axial Azimuthal coefficients vs r/R')
        plt.plot(r_R,Ct, label = r'$C_T$')
        plt.plot(r_R,Cq, label = r'$C_Q$')
        plt.grid()
        plt.xlabel('r/R')
        #plt.ylabel('Ct')
        plt.legend()
        plt.show()
        plt.savefig("Axial Azimuthal coefficients vs r_R.jpg", bbox_inches='tight')
        plt.close('all')
        
    @staticmethod
    def _plot_RotorThrustVSTipSpeedRatio(tip_speed_ratios,Rotor_Thrust):
        plt.figure(figsize=(12, 6))
        plt.title('Rotor Thrust vs Tip Speed Ratio')
        plt.plot(tip_speed_ratios,Rotor_Thrust)
       # plt.plot(r_R,Cq, label = r'$C_Q$')
        plt.grid()
        plt.xlabel('Tip Speed Ratio')
        plt.ylabel('Thrust (N)')
        # plt.legend()
        plt.show()
        plt.savefig("Rotor Thrust vs Tip Speed Ratio.jpg", bbox_inches='tight')
        plt.close('all')
    @staticmethod
    def _plot_RotorTorqueVSTipSpeedRatio(tip_speed_ratios,Rotor_Torque):
        plt.figure(figsize=(12, 6))
        plt.title('Rotor Torque vs Tip Speed Ratio')
        plt.plot(tip_speed_ratios,Rotor_Torque)
       # plt.plot(r_R,Cq, label = r'$C_Q$')
        plt.grid()
        plt.xlabel('Tip Speed Ratio')
        plt.ylabel('Torque (N.m)')
        # plt.legend()
        plt.show()
        plt.savefig("Rotor Torque vs Tip Speed Ratio.jpg", bbox_inches='tight')
        plt.close('all')   
        
        
        
        
        
        
        
        
        
        
        