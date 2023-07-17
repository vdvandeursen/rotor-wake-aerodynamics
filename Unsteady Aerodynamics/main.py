import math

from unsteady_panels import UnsteadyModel
import matplotlib.pyplot as plt


def plot_steady_cl_vs_alpha():
    # PLOTTING CL VS ALPHA FOR STEADY CASE
    lift_coefficients = []
    alpha_range = range(-15, 45, 10)
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
    panel_range = [1,2,5, 10, 20, 40,80]
    #step_range = [50, 100, 200]
    step_range=200
    res_dict = {}
    AoA = 10

    for npanels in panel_range:
        lift_coefficients = []
        print(f'running... {panel_range.index(npanels) + 1}/{len(panel_range)}')
        # for nsteps in step_range:
        model = UnsteadyModel()
        cl,cl_list, Tplot = model.run(
                alpha=AoA,
                reduced_freq=0.0,
                amplitude=0,
                dt=40,
                n_panels=npanels,
                n_steps=step_range,
                printing=False,
                plot_cl_t=False,
                plot_velocity_pressure=True
            )

        lift_coefficients=cl_list

        res_dict[npanels] = lift_coefficients

    plt.figure()

    for i in res_dict:
        plt.plot([i for i in range(step_range)], res_dict[i], label=f'{i} panels')
    plt.xlabel('Number of time steps', fontsize=12)
    plt.ylabel('Lift coefficient [-]', fontsize=12)
    plt.ylim([0.6,0.1+max([res_dict[i][-1] for i in res_dict])])
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f'./figures/cl_vs_npanels_vs_nsteps.png', dpi=600)

    plt.clf()
    plt.figure()
    for i in res_dict:
        plt.plot(Tplot, [x/(2*math.pi*AoA*math.pi/180.) for x in res_dict[i]], label=f'{i} panels')
    plt.plot(Tplot,[1-0.165*math.exp(-2*0.041*x)-0.335*math.exp(-2*0.3*x) for x in Tplot],'--', label="Wagner function")
    plt.xlabel('Distance travelled in chords', fontsize=12)
    plt.ylabel('Fraction of Steady Lift Coefficient [-]', fontsize=12)
    plt.ylim([0.4,1.05])
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f'./figures/cl_theory_vs_npanels_vs_nsteps.png', dpi=600)


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



def plot_cl_vs_k_vs_time():
    # PLOTTING DEPENDENCE ON N PANELS
    k_range = [0, 0.02, 0.05, 0.1]
    print()
    res_dict = {}
    o_dict = {}

    for k in k_range:
        print(f'running... {k_range.index(k) + 1}/{len(k_range)}')
        N=100
        if k!=0:
            i=2 #1 for pitch down, 2 for pitch up
            while N<101:
                N=round(i*math.pi*5/(k))+1
                print(N)
                i+=2
        model = UnsteadyModel()
        cl,cl_list,t_list, o_list = model.run(
                alpha=5,
                reduced_freq=k,
                amplitude=5,
                dt=5,
                n_panels=10,
                n_steps=N,
                printing=False,
                plot_cl_t=False,
                plot_velocity_pressure=True
            )

        res_dict[k] = cl_list
        o_dict[k] = o_list

    # plt.figure()
    #
    # for i in res_dict:
    #     plt.plot(t_list, res_dict[i], label=f'k={i}')
    # plt.xlabel('Distance travelled in chords', fontsize=12)
    # plt.ylabel('Lift coefficient [-]', fontsize=12)
    # plt.tight_layout()
    # plt.grid()
    # plt.legend()
    # plt.savefig(f'./figures/cl_vs_k_vs_time.png', dpi=600)
    # plt.clf()

    plt.figure()

    for i in res_dict:
        plt.plot([x*180/math.pi for x in o_dict[i][25:]], res_dict[i][25:], label=f'k={i}')
        plt.scatter(o_dict[i][-1] * 180 / math.pi, res_dict[i][-1])
    plt.plot([0,10],[0,2*math.pi*10*math.pi/180], "--", c="grey", label='Steady')
    plt.xlabel('AoA', fontsize=12)
    plt.ylabel('Lift coefficient [-]', fontsize=12)
    plt.ylim([-0.2,1.2])
    plt.xlim([-2,12])
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f'./figures/cl_vs_k_vs_aoa.png', dpi=600)

def plot_cl_vs_timestep():
    # PLOTTING DEPENDENCE ON N PANELS
    #panel_range = [1,2,5, 10, 20, 40,80]
    #step_range = [50, 100, 200]
    dts=[2,5,10,20,40]
    step_range=200
    res_dict = {}
    AoA = 10
    npanels = 20
    Tplot_dist={}

    for dt in dts:
        lift_coefficients = []
        print(f'running... {dts.index(dt) + 1}/{len(dts)}')
        # for nsteps in step_range:
        model = UnsteadyModel()
        cl,cl_list, Tplot = model.run(
                alpha=AoA,
                reduced_freq=0.0,
                amplitude=0,
                dt=dt,
                n_panels=npanels,
                n_steps=int(5*dt)+1,
                printing=False,
                plot_cl_t=False,
                plot_velocity_pressure=True
            )

        lift_coefficients=cl_list

        res_dict[dt] = lift_coefficients
        Tplot_dist[dt] = Tplot

    plt.figure()

    for i in res_dict:
        plt.plot(res_dict[i], label=f'dt=U/({i}C)')
    plt.xlabel('Number of time steps', fontsize=12)
    plt.ylabel('Lift coefficient [-]', fontsize=12)
    plt.ylim([0.6,0.1+max([res_dict[i][-1] for i in res_dict])])
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f'./figures/cl_vs_dt.png', dpi=600)

    plt.clf()
    plt.figure()
    for i in res_dict:
        plt.plot(Tplot_dist[i], [x/(2*math.pi*AoA*math.pi/180.) for x in res_dict[i]], label=f'dt=U/({i}C)')
    plt.plot(Tplot,[1-0.165*math.exp(-2*0.041*x)-0.335*math.exp(-2*0.3*x) for x in Tplot],'--', label="Wagner function",c='grey')
    plt.xlabel('Distance travelled in chords', fontsize=12)
    plt.ylabel('Fraction of Steady Lift Coefficient [-]', fontsize=12)
    plt.ylim([0.4,1.05])
    plt.tight_layout()
    plt.grid()
    plt.legend()
    plt.savefig(f'./figures/cl_theory_vs_dt.png', dpi=600)



if __name__ == '__main__':
    # Comment out the plots you do not want to run.

    #plot_steady_cl_vs_alpha()
    #plot_cl_vs_npanels_vs_nsteps()
    #plot_cl_vs_k_vs_ampl()
    #plot_cl_vs_timestep()
    plot_cl_vs_k_vs_time()

    print('break')
