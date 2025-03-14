import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from mpl_toolkits.mplot3d import Axes3D

def wigner_2d_dynamics(g1, kappa=0, gamma=0, detuning=0, n=0, N=15, times=[0, 5, 10, 15, 20, 25]):
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

    for t in times:
        idx = np.argmin(np.abs(t_list - t))  # Find closest time index
        rho_cavity = result.states[idx].ptrace(0)  # Trace out atomic state
        W = qt.wigner(rho_cavity, xvec, xvec)

        plt.figure(figsize=(8, 6))
        plt.contourf(xvec, xvec, W, 100, cmap='RdBu')
        plt.colorbar()
        plt.title(f'2D Wigner Function at t={t}')
        plt.xlabel('Re[$\\alpha$]')
        plt.ylabel('Im[$\\alpha$]')
        plt.show()


def wigner_3d_dynamics(g1, kappa=0, gamma=0, detuning=0, n=0, N=15, times=[0, 5, 10, 15, 20, 25]):
    """
    Computes and plots the 3D Wigner function for the Jaynes-Cummings model.

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

    for t in times:
        idx = np.argmin(np.abs(t_list - t))
        rho_cavity = result.states[idx].ptrace(0)
        W = qt.wigner(rho_cavity, xvec, xvec)

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        surf = ax.plot_surface(X, Y, W, rstride=1, cstride=1, cmap='RdBu', edgecolor='none')

        ax.set_xlabel('Re[$\\alpha$]')
        ax.set_ylabel('Im[$\\alpha$]')
        ax.set_zlabel('Wigner function')
        ax.set_title(f'3D Wigner Function at t={t}')
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()
