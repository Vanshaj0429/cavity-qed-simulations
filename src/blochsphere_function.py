import numpy as np
import qutip as qt
import matplotlib.pyplot as plt

def bloch_dynamics_snapshots(g1, kappa=0, gamma=0, detuning=0, n=0, N=15, timesteps=100, snapshot_times=[0, 5, 10, 15, 20, 25]):
    """
    Plots Bloch sphere snapshots at different time steps for Jaynes-Cummings model.

    Parameters:
    - g1: Coupling strength
    - kappa: Cavity decay rate
    - gamma: Atomic decay rate
    - detuning: Frequency detuning
    - n: Thermal photon number
    - N: Number of cavity states
    - timesteps: Number of time steps for visualization
    - snapshot_times: Specific times to visualize the state

    Returns:
    - Multiple Bloch spheres showing state evolution at different times.
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
    psi0 = qt.tensor(qt.basis(N, 0), qt.basis(2, 1))  # Cavity vacuum, atom excited
    t_list = np.linspace(0, 25, timesteps)
    result = qt.mesolve(H, psi0, t_list, c_ops, [])

    # Extract Bloch vector components
    bloch_vectors = np.zeros((3, len(t_list)))

    for i, state in enumerate(result.states):
        rho_atom = state.ptrace(1)  # Trace out cavity to get atomic density matrix
        bloch_vectors[:, i] = [
            qt.expect(qt.sigmax(), rho_atom),
            qt.expect(qt.sigmay(), rho_atom),
            qt.expect(qt.sigmaz(), rho_atom),
        ]

    # Plot Bloch spheres at selected time steps
    fig, axes = plt.subplots(1, len(snapshot_times), figsize=(18, 8), subplot_kw={'projection': '3d'})

    for i, t in enumerate(snapshot_times):
        idx = np.argmin(np.abs(t_list - t))  # Find closest time index
        bloch_sphere = qt.Bloch(fig=fig, axes=axes[i])
        
        # **Fix: Ensure correct shape**
        vector = np.array([[bloch_vectors[0, idx]], [bloch_vectors[1, idx]], [bloch_vectors[2, idx]]])
        
        bloch_sphere.add_vectors(vector.flatten())  # Correctly formatted 1D array
        bloch_sphere.add_points(bloch_vectors[:, :idx + 1])  # Correct shape (3, N)
        bloch_sphere.render()
        axes[i].set_title(f't = {t}')

    plt.show()
