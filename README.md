# CRYSTAL composite_methods
This repository contains a development release of the CRYSTAL code containing a new class of composite and low-cost electronic structure methods based on the r<sup>2</sup>SCAN functional, aimed at performing large-scale computations on modern HPC architectures for solid state calculations.

The corresponding r<sup>2</sup>SCAN-based methods, in their pure and hybrid HF/DFT variants, have been developed in combination with double-ζ (sol-def2-mSVP) and triple-ζ (pob-TZVP-rev2) quality basis sets specifically adapted to solids and augmented with the well-established gCP and D3(BJ) corrections.

The release includes a serial executable and the source code developed by us to compute the semi-empirical corrections exploited by these methods, together with some input examples. The executable contains a limitation w.r.t. number of atoms in the unit cell.

Further details can be found in the corresponding sections below.

## Executables
The CRYSTAL binary has been compiled for the Linux x86_64 architecture using Intel® oneAPI compilers.
- [Intel® oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html?operatingsystem=linux&linux-install-type=offline)
- [Intel® HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html?operatingsystem=linux&linux-install-type=offline)

## Running CRYSTAL with the new "sol-3c" and "pob-3c" composite methods
In the inputs folder, simply type the following command to test the new composite methods implemented in the CRYSTAL code, replacing input_name.d12 with the appropriate input file name.
```
./crystal < input_name.d12 >& output_name.out &
```

## License

This project is licensed under the MIT License
