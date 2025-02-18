import numpy as np
import matplotlib.pyplot as plt
import qutip as qt

def jaynes_cummings_dynamics(g1, kappa=0, gamma=0, detuning=0, n=0, non_linear=False):
    """
    Simulates the Jaynes-Cummings model using QuTiP.

    Parameters:
    - g1: Coupling strength.
    - kappa: Cavity decay rate.
    - gamma: Atomic decay rate.
    - detuning: Frequency detuning.
    - n: Average thermal photon number.
    - non_linear: If True, applies nonlinear corrections.

    Returns:
    - Plots of cavity and atomic occupation probabilities.
    - Power spectrum.
    """
    
    g = 2 * np.pi * g1
    omega_c_base = 2 * np.pi * 1  
    omega_a_base = 2 * np.pi * 1  

    gamma_phi = 0.001  # Dephasing rate

    if non_linear:
        omega_c = omega_c_base + detuning + 0.1 * g**2
        omega_a = omega_a_base + 0.1 * g**2
    else:
        omega_c = omega_c_base + detuning
        omega_a = omega_a_base

    # Operators
    N = 4  # Cavity states
    a = qt.tensor(qt.destroy(N), qt.qeye(2))
    sm = qt.tensor(qt.qeye(N), qt.destroy(2))

    c_ops = [
        np.sqrt(kappa * (1 + n)) * a, np.sqrt(kappa * n) * a.dag(),
        np.sqrt(gamma) * sm,
        np.sqrt(gamma_phi) * qt.tensor(qt.qeye(N), qt.sigmaz())
    ]

    H = omega_c * a.dag() * a + omega_a * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())

    # Initial state
    psi0 = qt.tensor(qt.basis(N, 0), qt.basis(2, 1))

    # Time evolution
    t_list = np.linspace(0, 25, 1000)
    result = qt.mesolve(H, psi0, t_list, c_ops, [a.dag() * a, sm.dag() * sm])

    # Plot Occupation Probabilities
    plt.figure(figsize=(8, 4))
    plt.plot(t_list, result.expect[0], label="Cavity")
    plt.plot(t_list, result.expect[1], label="Atom")
    plt.xlabel('Time (sec)')
    plt.ylabel('Occupation Probability')
    title = f'Jaynes-Cummings Dynamics: g={g1}, kappa={kappa}, gamma={gamma}, Detuning={detuning}, n={n}'
    plt.title(title + (", Non-linear" if non_linear else ""))
    plt.legend()
    plt.show()

    # Compute Correlation Function and Spectrum
    tlist = np.linspace(0, 100, 5000)
    corr = qt.correlation_2op_1t(H, None, tlist, c_ops, a.dag(), a)
    wlist1, spec1 = qt.spectrum_correlation_fft(tlist, corr)
    
    wlist2 = np.linspace(0.25, 1.75, 200) * 2 * np.pi
    spec2 = qt.spectrum(H, wlist2, c_ops, a.dag(), a)

    # Plot Power Spectrum
    plt.figure(figsize=(8,4))
    plt.plot(wlist1 / (2 * np.pi), spec1, 'b', lw=2, label='eseries method')
    plt.plot(wlist2 / (2 * np.pi), spec2, 'r--', lw=2, label='me+fft method')
    plt.xlabel('Frequency (meV)')
    plt.ylabel('Power Spectrum (arb. units)')
    plt.title(f'Vacuum Rabi Splitting: g={g1}, kappa={kappa}, gamma={gamma}, Detuning={detuning}, n={n}')
    plt.legend()
    plt.show()
