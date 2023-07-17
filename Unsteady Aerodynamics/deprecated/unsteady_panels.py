import math
import numpy as np
import matplotlib.pyplot as plt
from plotting_functions import plot_velocity_pressure_fields


class UnsteadyModel:
    def __init__(self):

        self.rho = 1.
        self.U_inf = 50.
        self.chord = 1.0
        # self.ALFA1 = 10.0
        self.vortices = None

    def run(self,
            alpha: float,
            reduced_freq: float,
            amplitude: float,
            n_panels: int,
            n_steps: int,
            dt: float,
            printing=True,
            plot_velocity_pressure=False,
            plot_cl_t=False):
        """
        Returns lift and circulation of vortex elements for a discretized flat plate in certain unsteady configuration.

        Parameters:
            :param alpha - angle of attack in deg.
            :param reduced_freq - reduced frequency of airfoil motion
            :param amplitude - amplitude of airfoil motion [deg]
            :param n_panels - number of panels to use to discretize airfoil.
            :param n_steps - number of timesteps to take.
            :param dt - Discretization timestep. Specified in Hz [1/s].
            :param printing - whether or not to print intermediate results
            :param plot_velocity_pressure - whether or not to plot the velocity and pressure field
            :param plot_cl_t - whether or not to plot the Cl vs time plot.

        """

        if printing:
            print(' SUDDEN ACCELERATION OF A FLAT PLATE\n' + 38 * '=' + '')
            print("  n      SX         CD         CL       GAMMAT       CLT")

        self.vortices = [[1.0] * 3 for _ in range(n_steps)]

        UW = [[0.0] * 2 for _ in range(n_steps)]  # Added UW initialization
        GAMAT1 = 0.0  # Added GAMAT1 initialization
        omega = reduced_freq * self.U_inf / self.chord  # oscillation frequency

        # self.DT = self.C / self.UT / 40.0
        dt = self.chord / self.U_inf / dt
        time = -dt

        CL_ratio = []
        T_plot = []

        alpha = np.radians(alpha)
        sine = math.sin(alpha)
        cosine = math.cos(alpha)

        DXW = 0.3 * self.U_inf * dt

        for n in range(1, n_steps + 1):
            time = time + dt
            amp = np.radians(amplitude)
            Omega = amp * math.sin(omega * time) + alpha
            Omega_dot = amp * omega * math.cos(omega * time)
            sine = math.sin(Omega)
            cosine = math.cos(Omega)

            # PATH OF ORIGIN (SX,SZ)
            SX = -self.U_inf * time
            SZ = 0.0

            # SHEDDING OF WAKE POINTS
            self.vortices[n - 1][0] = (self.chord + DXW) * cosine + SX
            self.vortices[n - 1][1] = -(self.chord + DXW) * sine + SZ
            # Create matrix containing a(ij)
            matrixA = np.ones([n_panels + 1, n_panels + 1])
            RHS = []

            for i in range(n_panels):  # collocation points
                # calculate collocation point
                x_cp = ((i + 0.75) * self.chord / n_panels) * cosine + SX
                z_cp = -((i + 0.75) * self.chord / n_panels) * sine + SZ
                for j in range(n_panels):  # panel vortex
                    x_j = ((j + 0.25) * self.chord / n_panels) * cosine + SX
                    z_j = -((j + 0.25) * self.chord / n_panels) * sine + SZ
                    U, W = self.VOR2D(x_cp, z_cp, x_j, z_j)
                    matrixA[i, j] = U * sine + W * cosine
                U, W = self.VOR2D(x_cp, z_cp, self.vortices[n - 1][0], self.vortices[n - 1][1])
                matrixA[i, n_panels] = U * sine + W * cosine

                U = math.cos(Omega) * self.U_inf
                W = math.sin(Omega) * self.U_inf + Omega_dot * x_cp
                Ww = 0
                if n != 1:
                    U1, W1 = self.DWASH(x_cp, z_cp, 1, n - 1)
                    Ww = W1 * math.cos(Omega) + U1 * math.sin(Omega)
                RHS.append(-(W + Ww))
            RHS2 = 0.0
            if n != 1:
                for I in range(1, n):
                    RHS2 = RHS2 - self.vortices[I - 1][2]
            RHS.append(RHS2)

            # CALCULATE MOMENTARY VORTEX STRENGTH OF WING AND WAKE VORTICES

            if n == 1:
                gammas1 = [0] * n_panels
            else:
                gammas1 = gammas.copy()

            gammas = np.matmul(np.linalg.inv(matrixA), RHS)
            self.vortices[n - 1][2] = gammas[-1]

            gammas = gammas[0:len(gammas) - 1]
            GAMMAT = sum(gammas)

            # WAKE ROLLUP
            if n >= 1:  # don't really understand the meaning of this conditional
                for I in range(1, n + 1):
                    utot = 0.0
                    wtot = 0.0
                    for j in range(n_panels):
                        x_j = ((j + 0.25) * self.chord / n_panels) * cosine + SX
                        z_j = -((j + 0.25) * self.chord / n_panels) * sine + SZ
                        U, W = self.VOR2D(self.vortices[I - 1][0], self.vortices[I - 1][1], x_j, z_j, gammas[j])
                        utot += U
                        wtot += W
                    U1, W1 = self.DWASH(self.vortices[I - 1][0], self.vortices[I - 1][1], 1, n)
                    U = utot + U1
                    W = wtot + W1
                    UW[I - 1][0] = self.vortices[I - 1][0] + U * dt  # Fixed UW assignment
                    UW[I - 1][1] = self.vortices[I - 1][1] + W * dt  # Fixed UW assignment

                for i in range(1, n + 1):
                    self.vortices[i - 1][0] = UW[i - 1][0]
                    self.vortices[i - 1][1] = UW[i - 1][1]

            # AERODYNAMIC LOADS
            if n == 1:
                GAMAT1 = 0.0

            dynamic_pressure = 0.5 * self.rho * self.U_inf ** 2
            DGAMDT = (GAMMAT - GAMAT1) / dt
            GAMAT1 = GAMMAT

            dgamsdt = [(gammas[j] - gammas1[j]) / dt for j in range(len(gammas))]
            dp = []
            for j in range(n_panels):
                dpj = self.rho * (self.U_inf * gammas[j] * n_panels / self.chord + sum(dgamsdt[0:j + 1]))
                dp.append(dpj)

            # CALCULATE WAKE INDUCED DOWNWASH
            XX1 = 0.75 * self.chord * cosine + SX
            ZZ1 = -0.75 * self.chord * sine + SZ
            U, W = self.DWASH(XX1, ZZ1, 1, n)
            WW = U * sine + W * cosine

            # L = self.rho * (self.U_inf * GAMMAT + DGAMDT * self.chord)
            # D = self.rho * (-WW * GAMMAT + DGAMDT * self.chord * sine)

            L = (sum([dp_ * self.chord / n_panels for dp_ in dp]))

            dd = []
            for i in range(n_panels):
                x_cp = ((i + 0.75) * self.chord / n_panels) * cosine + SX
                z_cp = -((i + 0.75) * self.chord / n_panels) * sine + SZ
                U, W = self.DWASH(x_cp, z_cp, 1, n)
                WW = U * sine + W * cosine
                ddj = self.rho * (-WW * gammas[i] + sum(dgamsdt) * self.chord / n_panels * sine)
                dd.append(ddj)

            D = sum(dd)

            CL = L / dynamic_pressure / self.chord
            CD = D / dynamic_pressure / self.chord

            CLT = CL / (2.0 * np.pi * sine)
            GAM1 = GAMMAT / (np.pi * self.U_inf * self.chord * sine)
            SX1 = SX - self.U_inf * dt

            CL_ratio.append(CL / ((Omega) * 2 * math.pi))
            T_plot.append(self.U_inf * time / self.chord)

            if printing:
                print(
                    f'{n:3}    {SX1:7.3f}    {CD:7.4f}    {CL:7.4f}    {GAM1:7.4f}    {CLT:7.4f} {(Omega) * 180 / math.pi:7.4f}    {CL / ((Omega) * 2 * math.pi):7.4f}  {self.U_inf * time / self.chord}')

        plot_name_suffix = f'{n_panels}panels_{n_steps}steps_{100*reduced_freq:.0f}k_{amplitude:.0f}ampl_{np.degrees(alpha):.0f}alpha'

        if plot_velocity_pressure:
            # Plot the velocity field: first determine vortex locations (plate + wake vortices)
            plate_vortex_elements = np.arange(n_panels)
            plate_vortex_x_locs = ((plate_vortex_elements + 0.25) * self.chord / n_panels) * cosine + SX
            plate_vortex_z_locs = -((plate_vortex_elements + 0.25) * self.chord / n_panels) * sine + SZ

            vortex_x_locs = np.concatenate([plate_vortex_x_locs, [vortex[0] for vortex in self.vortices]])
            vortex_z_locs = np.concatenate([plate_vortex_z_locs, [vortex[1] for vortex in self.vortices]])
            vortex_strengths = np.concatenate([gammas, [vortex[2] for vortex in self.vortices]])

            n = 3 * self.chord
            plot_velocity_pressure_fields(
                gammas=vortex_strengths,
                vortex_x_locs=vortex_x_locs,
                vortex_z_locs=vortex_z_locs,
                plate_start_loc=(SX, SZ),
                chord=self.chord,
                alpha=alpha,
                u_inf=self.U_inf,
                rho=self.rho,
                x_lim=(SX - n, SX + n),
                z_lim=(-n, n),
                n=(200, 100),
                plot_name_suffix=plot_name_suffix
            )

        if plot_cl_t:
            plt.figure()
            plt.plot(T_plot, CL_ratio)
            plt.ylabel('Lift coefficient [-]', fontsize=12)
            plt.xlabel('Nondimensional time [-]', fontsize=12)
            plt.ylim([0, 1.5])
            plt.tight_layout()
            plt.grid()
            plt.savefig(f'./figures/CL_T_plot_{plot_name_suffix}.png', dpi=600)
            # plt.show()

            plt.close()
        return CL

    def DWASH(self, X, Z, I1, I2):
        U = 0.0
        W = 0.0

        for I in range(I1, I2 + 1):
            U1, W1 = self.VOR2D(X, Z, self.vortices[I - 1][0], self.vortices[I - 1][1], self.vortices[I - 1][2])
            U = U + U1
            W = W + W1

        return U, W

    @staticmethod
    def VOR2D(X, Z, X1, Z1, GAMMA=1):
        U = 0.0
        W = 0.0
        RX = X - X1
        RZ = Z - Z1
        R = math.sqrt(RX ** 2 + RZ ** 2)

        if R < 0.001:
            return U, W

        V = 0.5 / np.pi * GAMMA / R
        U = V * (RZ / R)
        W = V * (-RX / R)

        return U, W
