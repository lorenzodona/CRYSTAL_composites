# composite_methods
This repository contains a development versions of the CRYSTAL code with a new class of composite electronic strucure methods based on the r<sup>2</sup>SCAN funtional aiming to perform large scale calculations on modern HPC architectures for ab initio solid-state materials simulations.

The first version, GCRYSTAL_1.0, directly integrates CUDA math libraries for accelerated linear algebra operations, building on the standard parallel version of the CRYSTAL code (PCRYSTAL). GCRYSTAL_1.0 has undergone extensive testing on accelerated HPC systems.

The second version, GCRYSTAL_2.0, is currently under development and involves a more extensive refactoring of the code. Its aim is to enhance the use of accelerated libraries by minimizing input/output operations and reducing frequent data transfers between the host and device. This version has not yet been thoroughly tested, and not all features are fully implemented.

Both versions are available as CRYSTAL DEMO executables (limited to 12 atoms), along with the source code for the newly developed code. Additional details can be found in the relative sections below.

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
