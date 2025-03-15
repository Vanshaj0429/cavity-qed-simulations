import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from mpl_toolkits.mplot3d import Axes3D

def wigner_2d_dynamics(g1, kappa=0, gamma=0, detuning=0, n=0, N=15, times=[0, 5, 10, 15, 20]):
    """
    Computes and plots the 2D Wigner function for the Jaynes-Cummings model.

    Parameters:
    - g1: Coupling strength
    - kappa: Cavity decay rate
    - gamma: Atomic decay rate
    - detuning: Frequency detuning
    - n: Thermal photon number
    - N: Number of cavity states
    - times: List of time steps for visualization
    """

    g = 2 * np.pi * g1
    wc = wa = 2 * np.pi * 1  # Base frequencies

    # Define operators
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.destroy(2))

    # Define Hamiltonian
    H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())

    # Define collapse operators
    c_ops = [
        np.sqrt(kappa * (1 + n)) * a,
        np.sqrt(kappa * n) * a.dag(),
        np.sqrt(gamma) * sm,
    ]

    # Solve the system
    psi0 = qt.tensor(qt.basis(N, 0), qt.basis(2, 1))
    t_list = np.linspace(0, 25, 100)
    result = qt.mesolve(H, psi0, t_list, c_ops, [])

    # Define phase space grid
    xvec = np.linspace(-5, 5, 200)

    # Create figure with subplots (VERTICAL layout)
    num_plots = len(times)
    fig, axes = plt.subplots(num_plots, 1, figsize=(8, 8 * num_plots), constrained_layout=False)

    if num_plots == 1:
        axes = [axes]  # Ensure axes is always a list for consistency

    # Generate Wigner functions and plot
    for i, t in enumerate(times):
        idx = np.argmin(np.abs(t_list - t))  # Find closest time index
        rho_cavity = result.states[idx].ptrace(0)  # Trace out atomic state
        W = qt.wigner(rho_cavity, xvec, xvec)

        ax = axes[i]
        contour = ax.contourf(xvec, xvec, W, 100, cmap='RdBu')

        ax.set_title(f't = {t}')
        ax.set_xlabel('Re[$\\alpha$]')
        ax.set_ylabel('Im[$\\alpha$]')

        # Add a colorbar for each subplot
        cbar = fig.colorbar(contour, ax=ax, fraction=0.05, pad=0.04)
        cbar.ax.tick_params(labelsize=8)  # Make colorbar labels smaller if needed

    # Adjust vertical spacing between plots
    plt.subplots_adjust(hspace=0.5)  # Increases vertical space between subplots
    plt.show()


def wigner_3d_dynamics(g1, kappa=0, gamma=0, detuning=0, n=0, N=15, times=[0, 5, 10, 15, 20]):
    """
    Computes and plots multiple 3D Wigner functions in a single figure.

    Parameters:
    - g1: Coupling strength
    - kappa: Cavity decay rate
    - gamma: Atomic decay rate
    - detuning: Frequency detuning
    - n: Thermal photon number
    - N: Number of cavity states
    - times: List of time steps for visualization
    """

    g = 2 * np.pi * g1
    wc = wa = 2 * np.pi * 1  # Base frequencies

    # Define operators
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.destroy(2))

    # Define Hamiltonian
    H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())

    # Define collapse operators
    c_ops = [
        np.sqrt(kappa * (1 + n)) * a,
        np.sqrt(kappa * n) * a.dag(),
        np.sqrt(gamma) * sm,
    ]

    # Solve the system
    psi0 = qt.tensor(qt.basis(N, 0), qt.basis(2, 1))
    t_list = np.linspace(0, 25, 100)
    result = qt.mesolve(H, psi0, t_list, c_ops, [])

    # Define phase space grid
    xvec = np.linspace(-5, 5, 50)
    X, Y = np.meshgrid(xvec, xvec)

    num_plots = len(times)

    # OPTION 1: Horizontal layout (1 row)
    #fig, axes = plt.subplots(1, num_plots, figsize=(num_plots * 5, 5), subplot_kw={'projection': '3d'})

    # OPTION 2: Vertical layout (1 column)
    fig, axes = plt.subplots(num_plots, 1, figsize=(10, num_plots * 10), subplot_kw={'projection': '3d'})

    if num_plots == 1:
        axes = [axes]  # Make sure axes is iterable when there's only one plot

    for i, t in enumerate(times):
        idx = np.argmin(np.abs(t_list - t))
        rho_cavity = result.states[idx].ptrace(0)
        W = qt.wigner(rho_cavity, xvec, xvec)

        surf = axes[i].plot_surface(X, Y, W, rstride=1, cstride=1, cmap='RdBu', edgecolor='none')

        axes[i].set_xlabel('Re[$\\alpha$]')
        axes[i].set_ylabel('Im[$\\alpha$]')
        axes[i].set_zlabel('Wigner function')
        axes[i].set_title(f't={t}')
        fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.tight_layout()
    plt.show()
