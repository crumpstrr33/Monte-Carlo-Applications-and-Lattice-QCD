**Note:** Looking over this code now, months after it was completed, I realize the myriad mistakes, inefficiencies, etc. that exist. If I were to write this code now, it would look much different and be high quality. Instead, I will keep it as is since this is my thesis work and fixing it now would be dishonest. I may in the future attempt again some of this code.
<br></br><br></br>

# Monte Carlo Applications and Lattice QCD
Contains all of the code used in the thesis along with the .pdf and .tex files and figures used for the thesis. Below are quick descriptions of the three main sections of the paper. Further explanations for the code can be found in the respective programs and for the physics of the code can be found in the paper (Honors Thesis.pdf) or in one of the individual sections (found in the Honors Thesis folder).

## Description
#### Ising Model
Calculated the effect of lattice size and temperature on the magnetization and susceptibility of a 2 dimensional Ising model. The lattice sizes investigated were
- N = 5
- N = 10
- N = 20
- N = 30
- N = 50

with temperatures from T/T<sub>c</sub> = 0.6 to 1.4 where T<sub>c</sub> is the critical temperature. The .mp4 file is a visual example of the Ising model that has a lattice size of 1000 produced by ising_animation.py.

#### Feynman Path Integral (on the Quantum Harmonic Oscillator)
After discretizing the Feynman path integral, it was applied to the ground state of the harmonic oscillator. Then it is applied to the correlation of an excited state where the excitation energy can be extracted. This was calculated by Monte Carlo where a large random sample size allows for a simple mean to approximate expectation values (in this case, the correlation function).

#### Gluonic Path Integral
Measured three different Wilson loops using a second- and third-order correction Wilson action as described in the paper. The Wilson loops are:
- AxA
- Ax2A
- Ax3A

where A is the lattice spacing. The Wilson loops are approximated via Monte Carlo analogously to the quantum harmonic oscialltor in the previous section.
