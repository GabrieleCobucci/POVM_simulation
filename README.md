# PVM-simulation of POVMs
Whereas standard quantum measurements are projective, the most general notion of a measurement is represented by positive operator-valued measures (POVMs). It is therefore natural to consider how accurately an experimenter with access only to projective measurements and classical processing can simulate POVMs. The most well-known class of non-projective measurements is called symmetric informationally complete (SIC). Such measurements are both ubiquitous in the broader scope of quantum information theory and known to be the most strongly non-projective measurements in qubit systems. Here, we show that beyond qubit systems, the SIC property is in general not associated with the most non-projective measurement. For this, we put forward a semidefinite programming criterion for detecting genuinely non-projective measurements. This method allows us to determine quantitative simulability thresholds for generic POVMs and to put forward a conjecture on which qutrit and ququart measurements that are most strongly non-projective.

# List of the branches:
- main: README.md
- Depolarising_noise: PVMsimulability.m, PVMsimulability_dual.m, most_non_proj.m
- Worst-noise_case: PVMsimulability_worst_noise.m, PVMsimulability_worst_noise_dual.m, most_non_proj_worst_noise.m
- Explicit_decomposition: check_decomposition.py, sic4_decomp.npy, sic5_decomp.npy

# List of the functions

- genpart.mlx : This function computes the rank-vectors associated to N-outcomes measurements in dimension d.

## Branch 1: Depolarising noise

- PVMsimulability.m : This function evaluates the critical visibility for PVM-simulability of a given POVM under depolarising noise via SDP;
- PVMsimulability_dual.m : This function evaluates the critical visibility for PVM-simulability of a given POVM under depolarising noise and witness operators to detect it;
- most_non_proj.m : This function finds the most non-projective N-outcome measurement in dimension d under depolarising noise.

## Branch 2: "Worst-noise" case

- PVMsimulability_worst_noise.m : This function evaluates the critical visibility for PVM-simulability of a given POVM under worst-case noise via SDP;
- PVMsimulability_worst_noise_dual.m : This function evaluates the critical visibility for PVM-simulability of a given POVM under worst-case noise and witness operators to detect it;
- most_non_proj_worst_noise.m : This function finds the most non-projective N-outcome measurement in dimension d under worst-case noise.

## Branch 3: Explicit_decomposition

- check_decomposition.py: This script checks whether a projective decomposition of a noisy SIC-POVM is valid;
- sic4_decomp.npy: Explicit decomposition for noisy SIC-POVM in d = 4;
- sic5_decomp.npy: Explicit decomposition for noisy SIC-POVM in d = 5.

# References

[1]: G. Cobucci, R. Brinster, S. Khandelwal, H. Kampermann, Dagmar Bruß, Nikolai Wyderka and Armin Tavakoli - Maximally non-projective measurements are not always symmetric informationally complete (https://arxiv.org/abs/2508.03652)
