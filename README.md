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

```bash
# Create the environment from the file
conda env create -f project.yml

# Activate the environment
conda activate [env_name_inside_yml]