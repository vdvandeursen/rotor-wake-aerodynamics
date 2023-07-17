import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def velocity_from_gammas(x, z, gammas, vortex_x_locs, vortex_z_locs):
    # V = gamma / (2pi distance)
    # Vx = V*cos theta with theta = atan(dz/dx)

    # thetas = np.arctan((z - vortex_z_locs)/(x - vortex_x_locs))
    thetas = np.arctan2((x - vortex_x_locs), (z - vortex_z_locs))
    distances = np.sqrt((x - vortex_x_locs) ** 2 + (z - vortex_z_locs) ** 2)

    # if min(distances) < 0.1:
    #     return np.nan, np.nan, np.nan

    Vx = sum(np.cos(thetas) * gammas / (2*np.pi*distances))
    Vz = - sum(np.sin(thetas) * gammas / (2*np.pi*distances))

    # if np.sqrt(Vx**2 + Vz**2) > 50:  # this happens when the point of interest is really close to a vortex.
    #     return np.nan, np.nan

    return Vx, Vz


def plot_velocity_pressure_fields(
        gammas:np.array,
        vortex_x_locs: np.array,
        vortex_z_locs: np.array,
        u_inf: float,
        rho: float,
        plate_start_loc: tuple,
        chord: float,
        alpha: float,
        x_lim: tuple,
        z_lim: tuple,
        n: tuple,
        plot_name_suffix=None):
    """Plots the velocity field around a superposition of vortex elements in a uniform flow.

    Parameters:
        :param gammas - array containing vortex strengths
        :param vortex_x_locs - vortex element chordwise coordinates
        :param vortex_z_locs - vortex element plungewise coordinates
        :param u_inf - free stream velocity
        :param rho - air density to calculate pressure field
        :param plate_start_loc - tuple with leading edge x,z coordinates
        :param chord - chord
        :param alpha - angle of attack
        :param x_lim - tuple wih start,end x-values of domain to plot
        :param z_lim - tuple with start,end z-values of domain to plot
        :param n - tuple with (nx, nz), the number of vectors to plot in x and z directions.
        :param plot_name_suffix - the name suffix.
    """
    if plot_name_suffix is None:
        plot_name_suffix = ""

    # First create the meshgrid of vectors to plot
    nx, nz = n
    x = np.linspace(x_lim[0], x_lim[1], nx)
    z = np.linspace(z_lim[0], z_lim[1], nz)
    xx, zz = np.meshgrid(x, z)

    # Evaluate velocities at all points in meshgrid using vectorized function.
    v_func = np.vectorize(velocity_from_gammas, excluded=[2, 3, 4])

    Vx_field, Vz_field = v_func(xx, zz, gammas, vortex_x_locs, vortex_z_locs)
    Vx_field = Vx_field + u_inf  # superimpose free stream velocity
    Vabs_field = np.sqrt((Vx_field**2 + Vz_field**2))
    P_field = 101325 - rho * (u_inf**2 - Vabs_field**2)

    # plotting the plate as a thick line
    plate_x = [plate_start_loc[0], plate_start_loc[0] + chord*np.cos(alpha)]
    plate_y = [plate_start_loc[1], plate_start_loc[1] - chord*np.sin(alpha)]

    # plot velocity field with streamlines
    fig, ax = plt.subplots(1, 1)
    plt.streamplot(xx, zz, Vx_field, Vz_field, color='grey', linewidth=0.5)  # plot streamlines
    contour = ax.contourf(xx, zz, Vabs_field, cmap='jet', levels=20)
    plt.plot(plate_x, plate_y, linewidth=3, color='k')
    cbar = fig.colorbar(contour)
    cbar.ax.set_ylabel('Velocity [m/s]')
    plt.tight_layout()
    plt.savefig(f'./figures/velocity/velocity_field_{plot_name_suffix}.png', dpi=600)
    plt.close()

    # plot pressure field
    fig, ax = plt.subplots(1, 1)
    contour = ax.contourf(xx, zz, P_field, levels=20)
    cbar = fig.colorbar(contour)
    cbar.ax.set_ylabel('Static pressure [Pa]')
    plt.plot(plate_x, plate_y, linewidth=3, color='r')
    plt.legend()
    # plt.imshow(P, origin='lower', interpolation='bilinear')
    plt.tight_layout()
    plt.savefig(f'./figures/pressure/pressure_field_{plot_name_suffix}.png', dpi=600)
    plt.close()

    return None
