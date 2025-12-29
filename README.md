# Divided Differences Simulation for Anharmonic Oscillators

This repository contains the C++ implementation of a novel simulation method based on **Divided Differences** and the **Permutation Matrix Representation (PMR)**.

This code was developed as part of the thesis **[Insert Thesis Title Here]** to simulate the dynamics of superconducting circuits driven by an **Anti-Symmetric (AS) pulse**. We demonstrate the method’s integration with the dispersive Hamiltonian model and analyze the resulting simulation dynamics.

Additionally, this project investigates the inherent **time-memory trade-off** of the method, examining the computational complexity and evaluating various optimization strategies to balance resource constraints.

## Table of Contents
- [Background & References](#background--references)
- [Results](#results)
- [Installation & Usage](#installation--usage)
- [Configuration](#configuration)
- [Future Work](#future-work)
- [Citation](#citation)

## Background & References
This work builds upon the foundational research by Asaf Diringer and Itay Hen.

1. **The Anti-Symmetric Pulse:** *Derivation of the pulse used to drive the system.* > [Insert Author, Title, Journal/Link (Year)]

2. **The Simulation Method:** *Derivation of the divided differences method implemented in this code.* > [Insert Author, Title, Journal/Link (Year)]

## Results
This repository provides the first successful implementation of the PMR simulation method for this specific Hamiltonian. We successfully replicated the expected theoretical dynamics for the AS pulse.

![Simulation Phase Space](path/to/phase_space_image.png)
*Fig 1: Phase space trajectory showing the conditional displacement achieved with the AS pulse.*

### Benchmarking
To validate the accuracy of the new method, we benchmarked our simulation results against the standard **QuTiP** solver. The comparison demonstrates strong agreement between the two methods.

![Comparison to QuTiP](path/to/comparison_image.png)
*Fig 2: Comparison of simulation accuracy against the QuTiP solver benchmark.*

In the thesis, we also derive the required **expansion order** ($q$) needed for convergence. This derivation allows for pre-calculation of computational costs before runtime.

> **Full Thesis:** [Link to Thesis PDF or University Repository]

## Installation & Usage

### Prerequisites
We recommend using the provided `project.yml` file to set up a consistent Conda environment for dependencies and compilation tools.

#### Create the environment from the file
conda env create -f project.yml

#### Activate the environment
conda activate [env_name_inside_yml]

### Compilation

We provide a shell script to handle compilation. To compile the simulation executable:
Bash

#### Make the script executable
chmod +x compile_multipiece_pulse.sh

#### Run the compilation script
./compile_multipiece_pulse.sh

Note: You can modify the target source file inside compile_multipiece_pulse.sh if needed.

#### Running the Simulation

Once the executable is compiled, you can run it directly. Ensure that your configuration file (params.json) is in the same directory (or the path specified in your code).

#### Bash

./[your_executable_name]

### Configuration

The simulation is fully configurable via the params.json file. This allows you to modify physical and computational parameters without recompiling the code.

#### Key Parameters:

* Physical Parameters: Define the Hamiltonian terms, pulse duration, and drive strength.

* Simulation Settings: Control the number of time steps, the expansion order (q), and the dimension of the Hilbert space.

* Output: The simulation outputs a JSON file containing the calculated transition amplitudes. The destination directory for this output is also defined within params.json.

## Future Work

This simulation technique utilizes the Permutation Matrix Representation (PMR) to solve quantum dynamics. While demonstrated here on superconducting circuits, the method is applicable to a broad range of quantum systems.

It is particularly effective for regimes characterized by:

* Short Evolution Times: Where the propagator does not require massive expansion orders.

* Weak Perturbations: Where the off-diagonal terms are small relative to the diagonal.

* Sparse Transition Graphs: Where the Hamiltonian connectivity (d) is low.

A key advantage of this framework is the analytic predictability of the computational cost. We provide a derivation that allows researchers to pre-estimate the required expansion order (q) based on the system parameters. This feature enables users to assess runtime feasibility and memory requirements before initiating costly high-performance computing tasks.