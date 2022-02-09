# slepc_gep_solver

A Sparse Generalized Eigenvalue Solver intended to be used with the [fem_2d library](https://github.com/jeremiah-corrado/fem_2d)

## Dependencies

* [CMake](https://cmake.org/)
* [MPIch](https://www.mpich.org/)
* [Petsc](https://petsc.org/release/overview/)
* [Slepc](https://slepc.upv.es/)

Linux is highly recomended, but it should be possible to install Slepc and Petsc on Windows or Mac if needed.

## Installation

1. Install Petsc. Follow the instructions [here](https://petsc.org/release/install/install_tutorial/)
2. Install Slepc (an extension of Petsc). Follow the instructions in the installation section of [this](https://petsc.org/release/install/install_tutorial/) pdf
3. Download this repository into your directory of choice
4. Navigate to your dicrectory in a terminal and run `make`. This should produce a binary file called **solve_gep**
5. Set the environment variable `GEP_SOLVE_DIR` to the directory where **solve_gep** is located. This will allow *fem_2d* to find the solver.

>The environment variable can be set perminantly on Linux by adding the line:
>GEP_SOLVE_DIR="the/actual/path/on/your/system/slepc_gep_solver/" 
>to your `.environment` file in `/etc/`


## Citation

Please include one or more of the following citations in any academic or commercial work based on this repository:

* ***metapaper***
* [Corrado, Jeremiah; Harmon, Jake; Notaros, Branislav (2021): A Refinement-by-Superposition Approach to Fully Anisotropic hp-Refinement for Improved Efficiency in CEM. TechRxiv. Preprint. https://doi.org/10.36227/techrxiv.16695163.v1](https://doi.org/10.36227/techrxiv.16695163.v1)
* [Harmon, Jake; Corrado, Jeremiah; Notaros, Branislav (2021): A Refinement-by-Superposition hp-Method for H(curl)- and H(div)-Conforming Discretizations. TechRxiv. Preprint. https://doi.org/10.36227/techrxiv.14807895.v1](https://doi.org/10.36227/techrxiv.14807895.v1)
