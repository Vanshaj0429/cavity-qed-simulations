import numpy as np
from matplotlib import pyplot as plt
import qutip

def power_spectrum_vs_coupling(N=15, wc=2*np.pi, wa=2*np.pi, kappa=0.7, gamma=0.2, n_th=0.25, g_values=None):
    """
    Computes and plots the power spectrum for different Jaynes-Cummings coupling strengths.

    Parameters:
    - N: Number of cavity Fock states
    - wc: Cavity frequency (in 2π units)
    - wa: Atom frequency (in 2π units)
    - kappa: Cavity dissipation rate
    - gamma: Atomic dissipation rate
    - n_th: Thermal photon number (average)
    - g_values: List or numpy array of coupling strengths to iterate over

    Returns:
    - Plot of power spectra for different coupling strengths
    """
    
    if g_values is None:
        g_values = np.linspace(0.01, 1.0, 6)  # Default coupling strengths

    tlist = np.linspace(0, 100, 5000)
    wlist2 = np.linspace(0.25, 1.75, 200) * 2 * np.pi

    plt.figure(figsize=(8, 4))

    for g_strength in g_values:
        g = g_strength * 2 * np.pi

        # Define the Jaynes-Cummings Hamiltonian
        a = qutip.tensor(qutip.destroy(N), qutip.qeye(2))
        sm = qutip.tensor(qutip.qeye(N), qutip.destroy(2))
        H = wc * a.dag() * a + wa * sm.dag() * sm + g * (a.dag() * sm + a * sm.dag())

        # Define collapse operators
        c_ops = [
            np.sqrt(kappa * (1 + n_th)) * a,
            np.sqrt(kappa * n_th) * a.dag(),
            np.sqrt(gamma) * sm,
        ]

        # Compute correlation function and spectrum
        corr = qutip.correlation_2op_1t(H, None, tlist, c_ops, a.dag(), a)
        wlist1, spec1 = qutip.spectrum_correlation_fft(tlist, corr)

        plt.plot(wlist1 / (2 * np.pi), spec1, lw=2, label=f'g={g_strength:.2f}')

    plt.legend()
    plt.xlabel('Frequency (meV)')
    plt.ylabel('Power Spectrum (arb. units)')
    plt.title(f'Vacuum Rabi Splitting for Different Coupling Strengths\n$\\kappa$={kappa}, $\\gamma$={gamma}')
    plt.xlim(-0.5, 2.5)
    plt.show()

# Run the function when executed as a script
if __name__ == "__main__":
    power_spectrum_vs_coupling()
