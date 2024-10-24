# composite_methods
This repository contains a development versions of the CRYSTAL code with a new class of composite and cost-effective electronic strucure methods based on the r<sup>2</sup>SCAN funtional aiming to perform large scale computations on modern HPC architectures for solid-state calculations.

The corresponding r<sup>2</sup>SCAN-based composite methods in their pure and hybrid HF/DFT variants have been developed in combination with double-ζ (sol-def2-mSVP) and triple-ζ (pob-TZVP-rev2) quality basis sets specifically adapted to solids and augmented with the well established gCP and D3(BJ) corrections.

This demo version of the CRYSTAL code is limited to 12 atoms and can be run in serial mode.
The repository also contains the source code with the semi-empirical corrections exploited by these methods along with some input cases.

Additional details can be found in the relative sections below.

## Prerequisites
Install Intel® oneAPI compilers for Linux x86_64 architecture
- [Intel® oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=linux&linux-install-type=offline)
- [ntel® HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=linux&linux-install-type=offline)

## Using CRYSTAL with the new "sol-3c" and "pob-3c" composite methods
To run the DEMO code the input file must be renamed “INPUT”, and the code must be executed in the directory of this file, as follows:
```
cp test_quartz.d12 INPUT
mpirun -np (number of process between 1 and 4) /path/to/GCRYSTAL_X.0_executables < INPUT >& test_quartz.out &
```

## License

This project is licensed under the MIT License
