from unsteady_panels import UnsteadyModel
import matplotlib.pyplot as plt


def plot_steady_cl_vs_alpha():
    # PLOTTING CL VS ALPHA FOR STEADY CASE
    lift_coefficients = []
    alpha_range = range(5, 45, 10)
    for alpha in alpha_range:
        print(f'running... {alpha_range.index(alpha) + 1}/{len(alpha_range)}')
        model = UnsteadyModel()
        cl = model.run(
            alpha=alpha,
            reduced_freq=0.0,
            amplitude=0,
            dt=40,
            n_panels=10,
            n_steps=200,
            printing=False,
            plot_cl_t=False,
            plot_velocity_pressure=True
        )

        lift_coefficients.append(cl)

    plt.figure()
    plt.plot(alpha_range, lift_coefficients, label='k=0')
    plt.xlabel('Angle of attack [deg]', fontsize=12)
    plt.ylabel('Lift coefficient [-]', fontsize=12)
    plt.tight_layout()
    plt.grid()
    plt.savefig(f'./figures/cl_vs_alpha_steady.png', dpi=600)


def plot_cl_vs_npanels_vs_nsteps():
    # PLOTTING DEPENDENCE ON N PANELS
    panel_range = [5, 10, 20, 50]
    step_range = [50, 100, 200]

    res_dict = {}

    for npanels in panel_range:
        lift_coefficients = []
        print(f'running... {panel_range.index(npanels) + 1}/{len(panel_range)}')
        for nsteps in step_range:
            model = UnsteadyModel()
            cl = model.run(
                alpha=10,
                reduced_freq=0.0,
                amplitude=0,
                dt=40,
                n_panels=npanels,
                n_steps=nsteps,
                printing=False,
                plot_cl_t=False,
                plot_velocity_pressure=True
            )

            lift_coefficients.append(cl)

        res_dict[npanels] = lift_coefficients

    plt.figure()

    for i in res_dict:
        plt.plot(step_range, res_dict[i], label=f'{i} panels')
    plt.xlabel('Number of time steps', fontsize=12)
    plt.ylabel('Lift coefficient [-]', fontsize=12)
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f'./figures/cl_vs_npanels_vs_nsteps.png', dpi=600)


def plot_cl_vs_k_vs_ampl():
    # PLOTTING DEPENDENCE ON N PANELS
    ampl_range = [0, 0.5, 1]
    k_range = [0, 0.02, 0.05, 0.1]

    res_dict = {}

    for ampl in ampl_range:
        lift_coefficients = []
        print(f'running... {ampl_range.index(ampl) + 1}/{len(ampl_range)}')
        for k in k_range:
            model = UnsteadyModel()
            cl = model.run(
                alpha=10,
                reduced_freq=k,
                amplitude=ampl,
                dt=40,
                n_panels=10,
                n_steps=100,
                printing=False,
                plot_cl_t=False,
                plot_velocity_pressure=True
            )

            lift_coefficients.append(cl)

        res_dict[ampl] = lift_coefficients

    plt.figure()

    for i in res_dict:
        plt.plot(k_range, res_dict[i], label=f'Ampl. {i} deg')
    plt.xlabel('Reduced frequency k [-]', fontsize=12)
    plt.ylabel('Lift coefficient [-]', fontsize=12)
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f'./figures/cl_vs_k_vs_ampl.png', dpi=600)



if __name__ == '__main__':
    # Comment out the plots you do not want to run.

    plot_steady_cl_vs_alpha()
    plot_cl_vs_npanels_vs_nsteps()
    plot_cl_vs_k_vs_ampl()

    print('break')
