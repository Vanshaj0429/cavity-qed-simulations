# Cavity QED Simulations

*A computational study of the Jaynes-Cummings Model with decoherence and dissipation.*

## Overview
This repository contains **numerical simulations** investigating the **Jaynes-Cummings Model (JCM)**, a fundamental framework in **Cavity Quantum Electrodynamics (QED)**. Our research explores **quantum coherence**, **Rabi oscillations**, and **Rabi splitting** under various coupling regimes, incorporating decoherence effects using the **Lindblad Master Equation**.

## Project Aims
This was a project I undertook at Cardiff University during my masters under supervison of **Dr Amy Morreau** and was structured to achieve the following aims:
- **Provide adaptable simulation tools** that enable users to use **custom JCH function** with adjustable **coupling strengths, decay rates, and other parameters** for flexible experimentation.
- **Analyze the impact of decoherence and dissipation** on quantum coherence, examining how different parameters such as **cavity loss, atomic decay, and dephasing** affect system dynamics.
- **Visualize fundamental quantum phenomena**, including **Vacuum Rabi Oscillations, Rabi Splitting, Wigner Functions, and Bloch Sphere dynamics**, to gain deeper insights into quantum behavior.
- **Facilitate research in quantum information science**, contributing to the understanding of **light-matter interactions** for applications in **quantum computing, quantum optics, and cavity QED systems**.


## Key Features
- **Jaynes-Cummings Model**: Two-level system interacting with a quantized electromagnetic field.
- **Lindblad Master Equation**: Models **cavity loss, atomic decay, and dephasing**.
- **Qutip-based Simulations**: Uses the **QuTiP** package for quantum dynamics.
- **Custom Hamiltonian Functions**: Predefined functions allow users to easily modify system parameters.
- **Custom Visualization tools**: Plots of **Vacuum Rabi Oscillations**, **Rabi Splitting**, **Bloch sphere evolution**, and **Wigner function representation** can be rapidly created just by choosing the parameters and putting them in the plotting functions.

## Installation
To run the simulations, ensure you have Python installed, then install the required dependencies:

```bash
pip install numpy matplotlib qutip
```

## Usage
After cloning the repository:

```bash
git clone https://github.com/Vanshaj0429/cavity-qed-simulations.git
cd cavity-qed-simulations
```

Run the Jupyter Notebook:
```bash
jupyter notebook
```

## Project Structure
```
📂 cavity-qed-simulations
│── 📜 README.md          # This file
│── 📜 cavity_qed.ipynb   # Jupyter Notebook (To be added)
│── 📜 requirements.txt   # Dependencies
│── 📜 plots/           # Contains generated plots
│── 📜 src/               # Python scripts and custom function for simulation
```
## Files & Functionality
This repository contains several scripts that provide different functionalities:

- **`jaynes_cummings.py`** – Simulates the Jaynes-Cummings dynamics and generates plots of **cavity and atomic occupation probabilities**, along with the **power spectrum** for different parameter values.
- **`power_spectrum_vs_coupling.py`** – Computes and plots the **power spectrum** for various **coupling strengths**, analyzing the effects of weak and strong coupling on quantum behavior.
- **`2d_3d_wigner_function.py`** – Generates **2D and 3D Wigner function visualizations**, offering insights into phase-space representations of quantum states over time.
- **`blochsphere_function.py`** – Plots **Bloch sphere snapshots** to illustrate the **state evolution of the atomic qubit**, showing coherence and decoherence effects at different time steps.

Each script is designed to be modular, allowing users to plug in different **decay rates, coupling strengths, and other system parameters** to analyze various configurations.

## Simulation Details
### 1. Theoretical Background
The **Jaynes-Cummings Model** describes the interaction between a two-level atom and a **single-mode cavity field**. The system is governed by the Hamiltonian:

$$
H = \hbar \omega_c a^\dagger a + \frac{1}{2} \hbar \omega_a \sigma_z + \hbar g (a^\dagger \sigma_- + a \sigma_+)
$$

where:
- $( \omega_c, \omega_a )$ are the photon and atom frequencies.
- $( a^\dagger, a )$ are the **creation** and **annihilation** operators.
- $ g $ is the **coupling strength**.

To introduce **decoherence**, we solve the **Lindblad Master Equation**:

$$
\frac{d\rho}{dt} = -\frac{i}{\hbar} [H, \rho] + \sum_j \left( L_j \rho L_j^\dagger - \frac{1}{2} \{ L_j^\dagger L_j, \rho \} \right)
$$

where $( L_j )$ are **Lindblad operators** modeling **cavity decay**, **atomic decay**, and **dephasing**.

### 2. Simulations & Results
- **Vacuum Rabi Oscillations**: Energy exchange between **atom** and **cavity field**.
- **Vacuum Rabi Splitting**: Spectral changes under **strong vs. weak coupling**.
- Plots and data are stored in the `/plots` directory.

## Contributing
We welcome contributions! If you'd like to improve the simulations or add features:
1. **Fork the repository**.
2. **Clone your fork** and create a new branch.
3. **Make changes and commit**.
4. **Submit a pull request**.

## References
- [QuTiP: Quantum Toolbox in Python](https://qutip.org/)
- Jaynes, E.T. & Cummings, F.W. (1963). *Comparison of quantum and semiclassical radiation theories*.
- [Lindblad Master Equation](https://qutip.org/docs/latest/guide/dynamics/dynamics-master.html)


