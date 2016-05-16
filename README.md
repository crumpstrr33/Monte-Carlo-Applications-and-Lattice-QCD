# Monte Carlo Applications and Lattice QCD
Contains all of the code used in the thesis along with the .pdf and .tex files and figures used for the thesis. Below are quick descriptions of the three main sections of the paper. Further explanations for the code can be found in the respective programs and for the topics/physics can be found in the paper (Honors Thesis.pdf) or in one of the individual sections (found in the Honors Thesis folder).

## Description
#### Ising Model
Calculated the effect of lattice size and temperature on the magnetization and susceptibility of a 2 dimensional Ising model. The .mp4 file is a visual example of the Ising model that has a lattice size of 1000 produced by ising_animation.py.

#### Feynman Path Integral (on the Quantum Harmonic Oscillator)
After discretizing the Feynman path integral, it was applied to the ground state of the harmonic oscillator. Then it is applied to the correlation of an excited state at different times where the excitation energy can be extracted.

#### Gluonic Path Integral
Measured three different Wilson loops using a second- and third-order correction Wilson action as described in the paper. The Wilson loops are:
- AxA
- Ax2A
- Ax3A

where A is the lattice spacing. The Wilson loops are approximated via Monte Carlo thus using a simple mean of the real traces of the loops.
